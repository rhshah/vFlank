"""Per-variant orchestration — the reusable use-case layer.

This is the presentation-free core of ``vflank small``: it turns a loaded MAF
DataFrame plus already-built sources into a stream of per-variant outcomes, and
collects that stream into a :class:`RunResult`. It contains **no** Typer/Rich/
``print``/file-writing — callers (the CLI, a future web service, a notebook)
render and persist. Diagnostics go through the ``vflank`` logger.

Layering: ``cli → pipeline → io → core``. ``pipeline`` may use ``io`` and
``core`` but must not import ``cli`` or any presentation library.

The design is a generator of discriminated outcomes (:class:`Processed` /
:class:`Skipped` / :class:`Duplicate`) plus an eager :func:`collect`. Progress is
the caller's concern — it falls out of iterating the generator — so the core
stays unaware of it. See ``docs/research/orchestration-extraction.md``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

import pandas as pd

from .core.consensus import BamConsensusSource, ConsensusFlankSource
from .core.flanks import FlankSource, ReferenceFlankSource
from .core.skips import categorize_skip
from .core.variant import Variant
from .errors import ConsensusError
from .io import emit_primer3 as primer3_io
from .io import fasta as fasta_io
from .io import maf as maf_io
from .io.maf import MafColumns
from .logging import get_logger

log = get_logger()


@dataclass(frozen=True, slots=True)
class Processed:
    """A variant that produced records."""

    variant: Variant
    records: list[str]                 # FASTA: raw + masked
    detail: dict                       # the per-variant report/summary row
    primer3: primer3_io.Primer3Record | None
    n_masked: int
    n_inserted: int
    used_consensus: bool
    flagged: bool
    truncation: str | None             # warning text, or None


@dataclass(frozen=True, slots=True)
class Skipped:
    """A row that could not become a variant or whose fetch failed."""

    message: str                       # human-readable, prefixed with the row index


@dataclass(frozen=True, slots=True)
class Duplicate:
    """A variant collapsed by dedup (same key already seen)."""

    variant: Variant


Outcome = Processed | Skipped | Duplicate


@dataclass
class RunResult:
    """Everything a caller needs to render a summary, write files, and report."""

    records: list[str] = field(default_factory=list)
    primer3: list[primer3_io.Primer3Record] = field(default_factory=list)
    rows: list[dict] = field(default_factory=list)            # per-variant detail
    skip_messages: list[str] = field(default_factory=list)
    skip_breakdown: Counter[str] = field(default_factory=Counter)  # categorised; .most_common()
    truncations: list[str] = field(default_factory=list)
    n_processed: int = 0
    n_skipped: int = 0
    n_duplicate: int = 0
    n_masked: int = 0
    n_consensus: int = 0
    n_flagged: int = 0
    n_inserted: int = 0


def iter_small(
    df: pd.DataFrame,
    *,
    cols: MafColumns,
    reference,
    gnomad,
    flank: int,
    af_threshold: float,
    dedup: bool,
    uppercase: bool,
    bam_resolver=None,
    policy=None,
    emit_primer3: bool = False,
    require_coverage: float = 0.0,
) -> Iterator[Outcome]:
    """Yield one :class:`Outcome` per MAF row.

    Owns the per-sample consensus cache and closes it on exit (``finally``), so
    the generator is safe to abandon early. ``reference`` and ``gnomad`` are
    built and owned by the caller (the caller closes them); this function only
    closes the consensus sources it lazily creates.

    With ``bam_resolver`` given (BAM mode) the flank source is the per-sample
    patient consensus, falling back to reference + gnomAD masking where no usable
    BAM exists for a sample — logged once per sample, never silently.
    """
    bam_mode = bam_resolver is not None
    ref_flank_source = ReferenceFlankSource(
        reference, gnomad, flank=flank, af_threshold=af_threshold
    )
    consensus_cache: dict[str, FlankSource] = {}
    bam_warned: set[str] = set()
    seen_variants: set[tuple] = set()

    def _source_for(variant: Variant) -> tuple[FlankSource, bool]:
        """Pick the flank source for a variant (consensus per-sample, or reference)."""
        if not bam_mode:
            return ref_flank_source, False
        sample = variant.sample
        if sample in consensus_cache:
            cached = consensus_cache[sample]
            return cached, not isinstance(cached, ReferenceFlankSource)
        bam_path = bam_resolver(sample)
        if bam_path is None:
            if sample not in bam_warned:
                log.warning("No BAM for sample %s — using reference + gnomAD masking", sample)
                bam_warned.add(sample)
            consensus_cache[sample] = ref_flank_source
            return ref_flank_source, False
        try:
            src: FlankSource = ConsensusFlankSource(
                reference, BamConsensusSource(bam_path, policy), gnomad,
                flank=flank, af_threshold=af_threshold,
            )
        except ConsensusError as exc:
            if sample not in bam_warned:
                log.warning("BAM unusable for %s (%s) — reference + gnomAD", sample, exc)
                bam_warned.add(sample)
            consensus_cache[sample] = ref_flank_source
            return ref_flank_source, False
        consensus_cache[sample] = src
        return src, True

    try:
        for row_idx, row in df.iterrows():
            variant, reason = maf_io.parse_variant_row(row, cols)
            if variant is None:  # reason is set by contract when variant is None
                yield Skipped(message=f"row {row_idx} {reason}")
                continue

            # Dedup by CHR_POS_REF_ALT (sample-independent reference/gnomAD
            # masking). With a BAM the consensus is patient-specific, so the key
            # also includes the sample -> one record per (variant, sample).
            if dedup:
                key: tuple = (variant.chrom, variant.start, variant.end,
                              variant.ref.upper(), variant.alt.upper())
                if bam_mode:
                    key = (*key, variant.sample)
                if key in seen_variants:
                    yield Duplicate(variant=variant)
                    continue
                seen_variants.add(key)

            source, used_consensus = _source_for(variant)
            try:
                fr = source.fetch(variant)
            except Exception as exc:  # noqa: BLE001  (any fetch failure -> reported skip)
                yield Skipped(
                    message=f"row {row_idx} {variant.gene} {variant.raw_chrom}:"
                            f"{variant.start} — fetch error: {exc}"
                )
                continue

            # Flag silently-truncated flanks. pysam returns a short string when a
            # window runs off a contig end rather than raising. A short left flank
            # near position 1 is expected; anything shorter than requested on the
            # right (or shorter than the position allows on the left) is a contig
            # boundary the user should know about — the record is still emitted.
            exp_left = min(flank, variant.start - 1)
            truncated = len(fr.left) < exp_left or len(fr.right) < flank
            truncation = (
                f"row {row_idx} {variant.gene} {variant.raw_chrom}:{variant.start} — "
                f"flank truncated near contig boundary "
                f"(left {len(fr.left)}/{exp_left}, right {len(fr.right)}/{flank})"
                if truncated else None
            )

            ref, alt = variant.ref, variant.alt
            if uppercase:
                fr = fr.upper()
                ref, alt = ref.upper(), alt.upper()

            sample_tag = variant.sample if bam_mode else None
            records = fasta_io.format_records(variant, fr, ref, alt, sample=sample_tag)
            primer3 = (
                primer3_io.small_variant_record(
                    fasta_io.record_id(variant, ref, alt, sample_tag),
                    fr.left, fr.right, fr.masked_left, fr.masked_right, ref,
                )
                if emit_primer3 else None
            )

            detail: dict = {
                "Sample": variant.sample[:32], "Gene": variant.gene,
                "Chrom": variant.raw_chrom, "Start": variant.start, "End": variant.end,
                "Ref": ref[:8], "Alt": alt[:8],
                "LeftLen": len(fr.left), "RightLen": len(fr.right),
                "NMasked": fr.n_masked, "NCorrected": fr.n_corrected, "Truncated": truncated,
            }
            flagged = False
            if bam_mode:
                source_label = "consensus" if used_consensus else "reference"
                covered_frac = (
                    fr.covered / fr.total
                    if (used_consensus and fr.total and fr.covered is not None)
                    else None
                )
                flagged = (
                    require_coverage > 0 and covered_frac is not None
                    and covered_frac < require_coverage
                )
                detail["Source"] = source_label
                detail["CoveredFrac"] = (
                    round(covered_frac, 3) if covered_frac is not None else ""
                )
                detail["NInserted"] = fr.inserted or 0
                detail["Flagged"] = flagged

            yield Processed(
                variant=variant, records=records, detail=detail, primer3=primer3,
                n_masked=fr.n_masked, n_inserted=fr.inserted or 0,
                used_consensus=used_consensus, flagged=flagged, truncation=truncation,
            )
    finally:
        for sample, src in consensus_cache.items():
            bam = getattr(src, "bam", None)
            if bam is not None:
                try:
                    bam.close()
                except Exception as exc:  # noqa: BLE001  (cleanup must not mask the run)
                    log.debug("Closing consensus BAM for %s failed: %s", sample, exc)


def collect(outcomes: Iterable[Outcome]) -> RunResult:
    """Accumulate an outcome stream into a :class:`RunResult` (records, rows,
    counts, and the categorised skip breakdown)."""
    r = RunResult()
    for o in outcomes:
        if isinstance(o, Processed):
            r.records.extend(o.records)
            if o.primer3 is not None:
                r.primer3.append(o.primer3)
            r.rows.append(o.detail)
            if o.truncation is not None:
                r.truncations.append(o.truncation)
            r.n_processed += 1
            r.n_masked += o.n_masked
            r.n_inserted += o.n_inserted
            r.n_consensus += int(o.used_consensus)
            r.n_flagged += int(o.flagged)
        elif isinstance(o, Skipped):
            r.skip_messages.append(o.message)
            r.n_skipped += 1
        else:  # Duplicate
            r.n_duplicate += 1
    # Categorise skips so a large, uniform skip set reads as one line per
    # category instead of a wall of identical messages.
    r.skip_breakdown = Counter(categorize_skip(m) for m in r.skip_messages)
    return r
