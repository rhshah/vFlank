"""Shared construction of the reference-sequence source for CLIs.

Both ``small`` and ``fusion`` choose a reference backend the same way — a local
indexed FASTA (default) or the UCSC API (no download) — so that logic lives in
one place, mirroring ``_masking.make_pop_source``.
"""

from __future__ import annotations

from pathlib import Path

from ..core.reference_api import ReferenceApiSource
from ..errors import VflankError
from ..io.reference import ReferenceFasta


def validate_ref_source(ref_source: str) -> None:
    if ref_source not in ("file", "api"):
        raise VflankError(f"--ref-source must be 'file' or 'api', got '{ref_source}'")


def make_reference_source(
    ref_source: str, ref_genome: Path | None, genome_build: str
) -> ReferenceFasta | ReferenceApiSource:
    """Build the reference source: a local FASTA (``file``) or the UCSC API (``api``).

    Raises ``VflankError`` if ``file`` is chosen without ``--ref-genome`` — never
    a silent fallback to one backend or the other.
    """
    if ref_source == "api":
        return ReferenceApiSource(genome_build)
    if ref_genome is None:
        raise VflankError(
            "--ref-genome is required with --ref-source file (the default). "
            "Pass an indexed FASTA, or use --ref-source api for the no-download UCSC backend."
        )
    return ReferenceFasta(ref_genome)
