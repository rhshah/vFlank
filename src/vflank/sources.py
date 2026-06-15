"""Construction of the reference and population-frequency sources from config.

These factories turn plain configuration (the ``--ref-source`` / ``--pop-source``
choices, a build, a directory) into the ports the pipeline consumes — a
presentation-free seam shared by the CLI and any other surface (the web service,
notebooks). They live in the library (not ``cli/``) precisely so
:func:`vflank.pipeline.run_small` can build sources without importing ``cli``.

Nothing here prints or parses arguments; failures raise :class:`VflankError`.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from .core.popfreq import GnomadStore
from .core.popfreq_api import GnomadApiSource
from .core.reference_api import ReferenceApiSource
from .errors import VflankError
from .io.reference import ReferenceFasta


def validate_ref_source(ref_source: str) -> None:
    if ref_source not in ("file", "api"):
        raise VflankError(f"--ref-source must be 'file' or 'api', got '{ref_source}'")


def validate_run_options(
    genome_build: str, ref_source: str, pop_source: str, pop_data: str,
    pop_vcf_dir: Path | None,
) -> None:
    """Validate the option set shared by ``small`` and ``fusion`` (CLI and the
    library ``run_*`` entrypoints both call this — one place, no duplication)."""
    if genome_build not in ("hg19", "hg38"):
        raise VflankError(f"--genome-build must be 'hg19' or 'hg38', got '{genome_build}'")
    validate_ref_source(ref_source)
    validate_pop_options(pop_source, pop_data)
    if pop_vcf_dir is not None and not pop_vcf_dir.is_dir():
        raise VflankError(f"--pop-vcf-dir is not a directory: {pop_vcf_dir}")


def make_reference_source(
    ref_source: str, ref_genome: Path | None, genome_build: str
) -> ReferenceFasta | ReferenceApiSource:
    """Build the reference source: a local FASTA (``file``) or the UCSC API (``api``).

    Raises :class:`VflankError` if ``file`` is chosen without a reference path —
    never a silent fallback to one backend or the other.
    """
    if ref_source == "api":
        return ReferenceApiSource(genome_build)
    if ref_genome is None:
        raise VflankError(
            "--ref-genome is required with --ref-source file (the default). "
            "Pass an indexed FASTA, or use --ref-source api for the no-download UCSC backend."
        )
    return ReferenceFasta(ref_genome)


def validate_pop_options(pop_source: str, pop_data: str) -> None:
    if pop_source not in ("vcf", "api"):
        raise VflankError(f"--pop-source must be 'vcf' or 'api', got '{pop_source}'")
    if pop_data not in ("genome", "exome", "both"):
        raise VflankError(f"--pop-data must be 'genome', 'exome', or 'both', got '{pop_data}'")


def make_pop_source(
    pop_source: str,
    pop_vcf_dir: Path | None,
    genome_build: str,
    pop_data: str,
    chroms: Iterable[str],
) -> GnomadStore | GnomadApiSource | None:
    """Build the masking source, or ``None`` if the VCF backend has no directory.

    For the VCF backend, runs ``preflight`` against ``chroms`` so a wholly-absent
    requested data kind fails fast (no silent genome-only fallback).
    """
    if pop_source == "api":
        return GnomadApiSource(genome_build, pop_data)
    if pop_vcf_dir is not None:
        store = GnomadStore(pop_vcf_dir, genome_build, pop_data)
        resolved = sorted(set(chroms))
        if resolved:
            store.preflight(resolved)
        return store
    return None
