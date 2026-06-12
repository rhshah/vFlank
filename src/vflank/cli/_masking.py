"""Shared construction of the population-frequency masking source for CLIs.

Both ``small`` and ``fusion`` select a masking backend the same way; this keeps
that logic in one place.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from ..core.popfreq import GnomadStore
from ..core.popfreq_api import GnomadApiSource
from ..errors import VflankError


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
    """Build the masking source, or None if VCF backend with no directory.

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
