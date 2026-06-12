"""Population allele-frequency masking source (gnomAD), local-VCF backend.

Masking can draw on gnomAD **genome** and/or **exome** data (``--pop-data``):
flanks often fall in non-coding regions where only genomes have data, while
exomes add power in coding regions. ``both`` masks the *union* (a position is
masked if it is a common SNP in either cohort).

The hot line-parsing kernel (:func:`parse_common_snp_positions`) is a pure
function over an iterable of raw VCF lines, so it is unit-testable without pysam
and is the natural seam to later swap for a Rust/noodles implementation.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from ..errors import PopFreqError
from ..logging import get_logger

log = get_logger()

# The two gnomAD data kinds we can mask against.
POP_DATA_KINDS = ("genome", "exome")

# gnomAD per-chromosome filename patterns, keyed by build then data kind.
# ``{chrom}`` is the bare chromosome (e.g. "1", "X"); the first existing match
# in each list wins. gnomAD ships per-chromosome files for both genomes and
# exomes on both builds.
GNOMAD_PATTERNS: dict[str, dict[str, list[str]]] = {
    "hg19": {
        "genome": [
            "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz",
            "gnomad.genomes.r2.1.2.sites.{chrom}.vcf.bgz",
        ],
        "exome": [
            "gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz",
        ],
    },
    "hg38": {
        "genome": [
            "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz",
            "gnomad.genomes.v4.0.sites.chr{chrom}.vcf.bgz",
        ],
        "exome": [
            "gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz",
            "gnomad.exomes.v4.0.sites.chr{chrom}.vcf.bgz",
        ],
    },
}


def kinds_for(pop_data: str) -> tuple[str, ...]:
    """Map a ``--pop-data`` value to the data kinds it consults.

    ``genome``/``exome`` -> that one; ``both`` -> both. Raises on anything else.
    """
    if pop_data == "both":
        return POP_DATA_KINDS
    if pop_data in POP_DATA_KINDS:
        return (pop_data,)
    raise ValueError(f"invalid pop_data {pop_data!r}; expected genome|exome|both")


def resolve_vcf_for_chrom(
    pop_vcf_dir: Path, bare_chrom: str, genome_build: str, data_kind: str
) -> Path | None:
    """Return the first existing gnomAD file for this chrom/build/kind, else None."""
    for pattern in GNOMAD_PATTERNS.get(genome_build, {}).get(data_kind, []):
        candidate = pop_vcf_dir / pattern.format(chrom=bare_chrom)
        if candidate.exists():
            return candidate
    return None


def example_filename(genome_build: str, data_kind: str) -> str:
    """A representative expected filename (chr1), for error/help messages."""
    patterns = GNOMAD_PATTERNS.get(genome_build, {}).get(data_kind, [])
    return patterns[0].format(chrom="1") if patterns else f"gnomad.{data_kind}s.*.vcf.bgz"


def build_chrom_vcf_map(
    pop_vcf_dir: Path, genome_build: str, chroms: list[str], data_kind: str
) -> dict[str, Path | None]:
    """Resolve the VCF (of one kind) for each chromosome (for coverage reports)."""
    return {c: resolve_vcf_for_chrom(pop_vcf_dir, c, genome_build, data_kind) for c in chroms}


def parse_common_snp_positions(rows: Iterable[str], af_threshold: float) -> list[int]:
    """Return 1-based positions of SNPs whose max AF/AF_grpmax >= threshold.

    Only single-base substitutions are considered (REF and every ALT length 1).
    ``rows`` is an iterable of raw tab-delimited VCF data lines. Works for both
    genome and exome gnomAD VCFs (identical INFO AF fields).
    """
    positions: list[int] = []
    for row in rows:
        fields = row.split("\t")
        if len(fields) < 8:
            continue

        ref = fields[3]
        if len(ref) != 1:
            continue
        alts = fields[4].split(",")
        if not all(len(a) == 1 and a not in (".", "*") for a in alts):
            continue

        af_values: list[float] = []
        for token in fields[7].split(";"):
            if token.startswith("AF=") or token.startswith("AF_grpmax="):
                for v in token.split("=", 1)[1].split(","):
                    try:
                        f = float(v)
                    except ValueError:
                        continue
                    if f == f:  # reject NaN (NaN != NaN)
                        af_values.append(f)

        if af_values and max(af_values) >= af_threshold:
            positions.append(int(fields[1]))
    return positions


class GnomadStore:
    """Lazy, per-(kind, chromosome) tabix cache over a directory of gnomAD VCFs.

    Honours ``--pop-data`` (genome / exome / both). For ``both``, queried
    positions are the union across the two cohorts. Tabix handles are opened on
    first use; each file's chr-notation is detected on open so queries use the
    file's own contig names.
    """

    def __init__(self, pop_vcf_dir: Path, genome_build: str, pop_data: str = "genome") -> None:
        self.dir = pop_vcf_dir
        self.build = genome_build
        self.kinds = kinds_for(pop_data)  # raises ValueError on bad input
        self._cache: dict[tuple[str, str], Any] = {}  # (kind, bare) -> TabixFile|None
        self._has_chr: dict[tuple[str, str], bool] = {}

    def preflight(self, chroms: Iterable[str]) -> None:
        """Fail fast if a requested data kind is wholly absent for these chroms.

        This is what prevents a silent genome-only fallback when ``--pop-data
        exome``/``both`` is requested without the exome files present. Per-chrom
        gaps (some chromosomes missing) remain warnings at query time.
        """
        chroms = list(chroms)
        missing = [
            kind for kind in self.kinds
            if not any(resolve_vcf_for_chrom(self.dir, c, self.build, kind) for c in chroms)
        ]
        if missing:
            raise PopFreqError(
                f"No {'/'.join(missing)} gnomAD VCF(s) found in {self.dir} for "
                f"build {self.build} (expected e.g. "
                f"{example_filename(self.build, missing[0])}). "
                f"Download them or choose a different --pop-data."
            )

    def _get(self, kind: str, bare: str):
        key = (kind, bare)
        if key in self._cache:
            return self._cache[key]

        path = resolve_vcf_for_chrom(self.dir, bare, self.build, kind)
        if path is None:
            log.warning(
                "No %s VCF for chr%s in %s — %s masking skipped for this chromosome",
                kind, bare, self.dir.name, kind,
            )
            self._cache[key] = None
            return None

        if not Path(str(path) + ".tbi").exists():
            log.warning(
                "TBI index missing: %s — %s masking skipped for chr%s (fix: tabix -p vcf %s)",
                path.name, kind, bare, path,
            )
            self._cache[key] = None
            return None

        try:
            import pysam

            tbx = pysam.TabixFile(str(path))
        except Exception as exc:  # noqa: BLE001 - surface as a warning, keep going
            log.warning("Could not open %s: %s", path.name, exc)
            self._cache[key] = None
            return None

        self._cache[key] = tbx
        self._has_chr[key] = any(c.startswith("chr") for c in tbx.contigs)
        log.debug(
            "Opened %s %s (contigs: %s)",
            kind, path.name, "chr-prefixed" if self._has_chr[key] else "bare",
        )
        return tbx

    def get_positions(
        self, bare: str, start_0based: int, end_0based: int, af_threshold: float
    ) -> list[int]:
        """1-based positions of common SNPs in ``[start, end)`` — union over kinds."""
        positions: set[int] = set()
        for kind in self.kinds:
            tbx = self._get(kind, bare)
            if tbx is None:
                continue
            contig = f"chr{bare}" if self._has_chr.get((kind, bare)) else bare
            try:
                rows = tbx.fetch(contig, start_0based, end_0based)
            except (ValueError, KeyError):
                # Contig absent in this VCF (e.g. MT). Expected often enough to be
                # DEBUG, not WARNING — but never wholly silent.
                log.debug(
                    "Contig %r not in %s VCF for chr%s — no masking for this region",
                    contig, kind, bare,
                )
                continue
            positions.update(parse_common_snp_positions(rows, af_threshold))
        return sorted(positions)

    def close(self) -> None:
        for tbx in self._cache.values():
            if tbx is not None:
                tbx.close()
