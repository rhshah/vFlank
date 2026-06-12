"""Reference FASTA access with chr-notation detection and build fingerprinting.

The build-mismatch guard addresses the scariest silent failure in this domain:
running hg19 coordinates against an hg38 FASTA (or vice versa) returns the wrong
sequence with no error. We fingerprint by the length of chromosome 1.
"""

from __future__ import annotations

from pathlib import Path

from ..core.chrom import chrom_for_contigs, contigs_have_chr
from ..errors import ReferenceError

# chr1 length uniquely distinguishes the two common human builds.
_CHR1_LENGTH_TO_BUILD: dict[int, str] = {
    249_250_621: "hg19",  # GRCh37
    248_956_422: "hg38",  # GRCh38
}


class ReferenceFasta:
    """Thin wrapper over ``pysam.FastaFile`` keyed by bare chromosome."""

    def __init__(self, path: Path | str):
        self.path = Path(path)
        if not self.path.exists():
            raise ReferenceError(f"Reference FASTA not found: {self.path}")
        fai = Path(str(self.path) + ".fai")
        if not fai.exists():
            raise ReferenceError(
                f"FASTA index not found: {fai}  (fix: samtools faidx {self.path})"
            )
        try:
            import pysam

            self.fa = pysam.FastaFile(str(self.path))
        except Exception as exc:  # noqa: BLE001
            raise ReferenceError(f"Could not open FASTA: {exc}") from exc

        self._refs = set(self.fa.references)
        self.has_chr = contigs_have_chr(self.fa.references)

    def contig(self, bare: str) -> str:
        """Resolve a bare chromosome to this FASTA's actual contig name.

        Auto-detection of ``chr`` prefixing is best-effort (it probes chr1-5).
        If the detected form is absent but the other form is present — which can
        happen with unusual or single-contig references — fall back to it rather
        than letting pysam raise a confusing KeyError for every variant.
        """
        primary = chrom_for_contigs(bare, self.has_chr)
        if primary in self._refs:
            return primary
        alternate = bare if self.has_chr else f"chr{bare}"
        return alternate if alternate in self._refs else primary

    def fetch(self, bare: str, start_0based: int, end_0based: int) -> str:
        return self.fa.fetch(self.contig(bare), start_0based, end_0based)

    def detect_build(self) -> str | None:
        """Infer 'hg19'/'hg38' from chr1 length, or None if undeterminable."""
        for name in ("chr1", "1"):
            if name in self.fa.references:
                return _CHR1_LENGTH_TO_BUILD.get(self.fa.get_reference_length(name))
        return None

    def check_build(self, declared: str) -> str | None:
        """Return a warning string if the declared build disagrees with the FASTA."""
        detected = self.detect_build()
        if detected is not None and detected != declared:
            return (
                f"Declared --genome-build {declared} but the FASTA looks like "
                f"{detected} (by chr1 length). Coordinates may map to the wrong "
                f"sequence. Double-check the reference."
            )
        return None

    def close(self) -> None:
        self.fa.close()
