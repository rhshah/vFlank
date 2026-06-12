"""Population allele-frequency masking source (gnomAD).

The hot line-parsing kernel (:func:`parse_common_snp_positions`) is a pure
function over an iterable of raw VCF lines, so it is unit-testable without pysam
and is the natural seam to later swap for a Rust/noodles implementation.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from ..logging import get_logger

log = get_logger()

# gnomAD filename patterns, keyed by genome build. ``{chrom}`` is the bare
# chromosome (e.g. "1", "X"). Patterns without ``{chrom}`` are single-file
# resources covering all chromosomes (legacy fallback).
GNOMAD_PATTERNS: dict[str, list[str]] = {
    "hg38": [
        "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz",
        "gnomad.genomes.v4.0.sites.chr{chrom}.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr{chrom}.vcf.bgz",
        "gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz",
    ],
    "hg19": [
        "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz",
        "gnomad.genomes.r2.1.2.sites.{chrom}.vcf.bgz",
        "ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz",
    ],
}


def resolve_vcf_for_chrom(
    pop_vcf_dir: Path, bare_chrom: str, genome_build: str
) -> Path | None:
    """Return the first gnomAD file matching a known pattern for this chrom."""
    for pattern in GNOMAD_PATTERNS.get(genome_build, []):
        if "{chrom}" not in pattern:
            candidate = pop_vcf_dir / pattern
        else:
            candidate = pop_vcf_dir / pattern.format(chrom=bare_chrom)
        if candidate.exists():
            return candidate
    return None


def build_chrom_vcf_map(
    pop_vcf_dir: Path, genome_build: str, chroms: list[str]
) -> dict[str, Path | None]:
    return {c: resolve_vcf_for_chrom(pop_vcf_dir, c, genome_build) for c in chroms}


def parse_common_snp_positions(
    rows: Iterable[str], af_threshold: float
) -> list[int]:
    """Return 1-based positions of SNPs whose max AF/AF_grpmax >= threshold.

    Only single-base substitutions are considered (REF and every ALT length 1).
    ``rows`` is an iterable of raw tab-delimited VCF data lines.
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
    """Lazy, per-chromosome tabix cache over a directory of gnomAD VCFs.

    Tabix handles are opened on first use and reused; each file's chr-notation
    is detected on open so queries always use the file's own contig names.
    """

    def __init__(self, pop_vcf_dir: Path, genome_build: str) -> None:
        self.dir = pop_vcf_dir
        self.build = genome_build
        self._cache: dict[str, Any] = {}  # bare chrom -> pysam.TabixFile | None (untyped)
        self._has_chr: dict[str, bool] = {}

    def _get(self, bare: str):
        if bare in self._cache:
            return self._cache[bare]

        path = resolve_vcf_for_chrom(self.dir, bare, self.build)
        if path is None:
            log.warning("No VCF for chr%s in %s — masking skipped here", bare, self.dir.name)
            self._cache[bare] = None
            return None

        if not Path(str(path) + ".tbi").exists():
            log.warning(
                "TBI index missing: %s — masking skipped for chr%s "
                "(fix: tabix -p vcf %s)",
                path.name, bare, path,
            )
            self._cache[bare] = None
            return None

        try:
            import pysam

            tbx = pysam.TabixFile(str(path))
        except Exception as exc:  # noqa: BLE001 - surface as a warning, keep going
            log.warning("Could not open %s: %s", path.name, exc)
            self._cache[bare] = None
            return None

        self._cache[bare] = tbx
        self._has_chr[bare] = any(c.startswith("chr") for c in tbx.contigs)
        log.debug(
            "Opened %s (contigs: %s)",
            path.name, "chr-prefixed" if self._has_chr[bare] else "bare",
        )
        return tbx

    def get_positions(
        self, bare: str, start_0based: int, end_0based: int, af_threshold: float
    ) -> list[int]:
        """1-based positions of common SNPs in ``[start, end)`` for this chrom."""
        tbx = self._get(bare)
        if tbx is None:
            return []
        contig = f"chr{bare}" if self._has_chr.get(bare) else bare
        try:
            rows = tbx.fetch(contig, start_0based, end_0based)
        except (ValueError, KeyError):
            # Contig absent in this VCF (e.g. MT in some builds). Expected often
            # enough that it is DEBUG, not WARNING — but never wholly silent.
            log.debug("Contig %r not in VCF for chr%s — no masking for this region", contig, bare)
            return []
        return parse_common_snp_positions(rows, af_threshold)

    def close(self) -> None:
        for tbx in self._cache.values():
            if tbx is not None:
                tbx.close()
