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
from pathlib import Path

import pandas as pd

from .core.chrom import normalise_chrom
from .core.consensus import BamConsensusSource, ConsensusFlankSource
from .core.flanks import FlankSource, ReferenceFlankSource
from .core.fusion import build_junction
from .core.skips import categorize_skip
from .core.variant import Variant
from .errors import ConsensusError, VflankError
from .io import breakpoints as bp_io
from .io import emit_primer3 as primer3_io
from .io import fasta as fasta_io
from .io import maf as maf_io
from .io import vcf as vcf_io
from .io.breakpoints import SvColumns, SvInput
from .io.maf import MAF_CHR, MafColumns, MafInput
from .logging import get_logger
from .sources import make_pop_source, make_reference_source, validate_run_options

log = get_logger()

__all__ = [
    # batteries-included entrypoints (build sources -> run -> RunResult)
    "run_small", "run_fusion",
    # streaming primitives (caller drives progress / accumulation)
    "iter_small", "iter_fusion", "collect",
    # data types
    "RunResult", "Processed", "Skipped", "Duplicate", "Outcome",
]


@dataclass(frozen=True, slots=True)
class Processed:
    """A variant or fusion that produced records (raw + masked)."""

    records: list[str]                 # FASTA: raw + masked
    detail: dict                       # the per-item report/summary row
    primer3: primer3_io.Primer3Record | None
    n_masked: int
    used_consensus: bool
    n_inserted: int = 0                # small/BAM only
    flagged: bool = False              # small/--require-coverage only
    truncation: str | None = None      # small only: per-variant truncation message


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
    ref_api_requests: int | None = None   # UCSC reference API calls (--ref-source api)
    api_requests: int | None = None       # gnomAD API calls (--pop-source api)


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
                records=records, detail=detail, primer3=primer3,
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


def iter_fusion(
    df: pd.DataFrame,
    *,
    cols: SvColumns,
    reference,
    gnomad,
    flank: int,
    af_threshold: float,
    bam_resolver=None,
    policy=None,
    emit_primer3: bool = False,
) -> Iterator[Outcome]:
    """Yield one :class:`Outcome` per breakpoint row (:class:`Processed` or
    :class:`Skipped` — fusions are not deduplicated).

    Owns the per-sample BAM cache and closes it on exit (``finally``). Where a
    sample has no usable BAM the masking falls back to gnomAD only — logged once
    per sample, never silently.
    """
    bam_mode = bam_resolver is not None
    bam_cache: dict[str, BamConsensusSource | None] = {}
    bam_warned: set[str] = set()

    def _bam_for(sample: str) -> BamConsensusSource | None:
        if not bam_mode or not sample:
            return None
        if sample in bam_cache:
            return bam_cache[sample]
        path = bam_resolver(sample)
        src: BamConsensusSource | None = None
        if path is not None:
            try:
                src = BamConsensusSource(path, policy)
            except ConsensusError as exc:
                if sample not in bam_warned:
                    log.warning("BAM unusable for %s (%s) — gnomAD masking only", sample, exc)
                    bam_warned.add(sample)
        elif sample not in bam_warned:
            log.warning("No BAM for sample %s — gnomAD masking only", sample)
            bam_warned.add(sample)
        bam_cache[sample] = src
        return src

    try:
        for row_idx, row in df.iterrows():
            fusion, reason = bp_io.parse_fusion_row(row, cols)
            if fusion is None:  # reason is set by contract when fusion is None
                yield Skipped(message=f"row {row_idx} — {reason}")
                continue
            bam_source = _bam_for(fusion.sample)
            try:
                jr = build_junction(
                    reference, fusion, flank, gnomad=gnomad,
                    af_threshold=af_threshold, bam_source=bam_source,
                )
            except Exception as exc:  # noqa: BLE001  (any junction failure -> reported skip)
                yield Skipped(message=f"row {row_idx} {fusion.name} — junction error: {exc}")
                continue

            label = fasta_io.safe_header(fusion.name or "fusion")
            bp = (f"{fusion.bp1.chrom}_{fusion.bp1.pos}_{fusion.bp1.strand}__"
                  f"{fusion.bp2.chrom}_{fusion.bp2.pos}_{fusion.bp2.strand}")
            prefix = f"{fasta_io.safe_header(fusion.sample)}__" if fusion.sample else ""
            header = f"{prefix}{label}__{bp}__j{jr.junction_index}"
            records = [
                f">{header}\n{jr.sequence}\n",
                f">Masked__{header}\n{jr.masked_sequence}\n",
            ]
            primer3 = (
                primer3_io.junction_record(
                    header, jr.sequence, jr.masked_sequence, jr.junction_index
                )
                if emit_primer3 else None
            )
            detail = {
                "Name": fusion.name or ".", "BP1": f"{fusion.bp1.chrom}:{fusion.bp1.pos}",
                "BP2": f"{fusion.bp2.chrom}:{fusion.bp2.pos}",
                "Len": len(jr.sequence), "Junction": jr.junction_index,
                "N": jr.n_masked, "Trunc": len(jr.sequence) < 2 * flank,
            }
            yield Processed(
                records=records, detail=detail, primer3=primer3,
                n_masked=jr.n_masked, used_consensus=bam_source is not None,
            )
    finally:
        for sample, src in bam_cache.items():
            if src is not None:
                try:
                    src.close()
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


# --- Batteries-included entrypoints ---------------------------------------
# build sources from config -> run -> close -> RunResult. These are what a web
# service / notebook calls; the CLI keeps its own streaming flow (for the
# progress bar and per-source status lines) but shares iter_small/iter_fusion,
# the source factories, and validate_run_options with these.


def run_small(
    maf: MafInput,
    *,
    genome_build: str,
    cols: MafColumns | None = None,
    ref_source: str = "file",
    ref_genome: Path | None = None,
    pop_source: str = "vcf",
    pop_vcf_dir: Path | None = None,
    pop_data: str = "genome",
    flank: int = 200,
    af_threshold: float = 0.001,
    dedup: bool = True,
    uppercase: bool = True,
    sample_filter: set[str] | None = None,
    emit_primer3: bool = False,
) -> RunResult:
    """Run the small-variant pipeline end to end and return a :class:`RunResult`.

    Builds the reference and (optional) gnomAD sources from the given options,
    loads the input, runs the per-variant orchestration, closes the sources, and
    returns records + per-variant rows + skips + counts. The input is a MAF or a
    VCF/BCF: a VCF *path* is auto-detected by extension and read sites-only
    (``cols`` overrides do not apply); an open buffer always reads as MAF.
    Presentation-free: raises :class:`VflankError` on bad options; never prints.
    BAM consensus is a CLI-only path for now and is not exposed here.
    """
    cols = cols or MafColumns()
    validate_run_options(genome_build, ref_source, pop_source, pop_data, pop_vcf_dir)

    reference = make_reference_source(ref_source, ref_genome, genome_build)
    gnomad = None
    try:
        build_warn = reference.check_build(genome_build)
        if build_warn:
            log.warning("%s", build_warn)

        if vcf_io.is_vcf_path(maf):
            df = vcf_io.load_vcf(maf)
            if sample_filter is not None:
                log.warning("--samples ignored for VCF input (sites-only).")
                sample_filter = None
        else:
            df = maf_io.load_maf(maf, cols)
        if sample_filter is not None:
            if cols.sample not in df.columns:
                raise VflankError(f"Sample column '{cols.sample}' not found; cannot filter.")
            df = df[df[cols.sample].astype(str).isin(sample_filter)].copy()
            if df.empty:
                raise VflankError("No variants remain after sample filtering.")

        maf_chroms = {
            b for b, _err in (normalise_chrom(c) for c in df[MAF_CHR].dropna().unique()) if b
        }
        gnomad = make_pop_source(pop_source, pop_vcf_dir, genome_build, pop_data, maf_chroms)

        result = collect(iter_small(
            df, cols=cols, reference=reference, gnomad=gnomad, flank=flank,
            af_threshold=af_threshold, dedup=dedup, uppercase=uppercase,
            emit_primer3=emit_primer3,
        ))
        result.ref_api_requests = (
            getattr(reference, "request_count", None) if ref_source == "api" else None
        )
        result.api_requests = getattr(gnomad, "request_count", None) if gnomad is not None else None
        return result
    finally:
        reference.close()
        if gnomad is not None:
            gnomad.close()


def run_fusion(
    sv_table: SvInput,
    *,
    genome_build: str,
    cols: SvColumns | None = None,
    ref_source: str = "file",
    ref_genome: Path | None = None,
    pop_source: str = "vcf",
    pop_vcf_dir: Path | None = None,
    pop_data: str = "genome",
    flank: int = 200,
    af_threshold: float = 0.001,
    emit_primer3: bool = False,
) -> RunResult:
    """Run the fusion pipeline end to end and return a :class:`RunResult`.

    Mirrors :func:`run_small`: builds sources, loads the breakpoint table (path
    or buffer), builds junctions, closes sources. Presentation-free.
    """
    cols = cols or SvColumns()
    validate_run_options(genome_build, ref_source, pop_source, pop_data, pop_vcf_dir)

    reference = make_reference_source(ref_source, ref_genome, genome_build)
    gnomad = None
    try:
        build_warn = reference.check_build(genome_build)
        if build_warn:
            log.warning("%s", build_warn)

        df = bp_io.load_sv_table(sv_table, cols)
        bp_chroms = {
            b
            for col in (cols.chr1, cols.chr2)
            for b, _err in (normalise_chrom(v) for v in df[col].dropna().unique())
            if b
        }
        gnomad = make_pop_source(pop_source, pop_vcf_dir, genome_build, pop_data, bp_chroms)

        result = collect(iter_fusion(
            df, cols=cols, reference=reference, gnomad=gnomad, flank=flank,
            af_threshold=af_threshold, emit_primer3=emit_primer3,
        ))
        result.ref_api_requests = (
            getattr(reference, "request_count", None) if ref_source == "api" else None
        )
        result.api_requests = getattr(gnomad, "request_count", None) if gnomad is not None else None
        return result
    finally:
        reference.close()
        if gnomad is not None:
            gnomad.close()
