"""Chromosome notation detection and normalisation.

All functions here are pure (no I/O) so they are trivially unit-testable. The
canonical internal form is the *bare* chromosome string (e.g. ``"7"``, ``"X"``,
``"MT"``); callers convert to ``chr``-prefixed form for a given FASTA/VCF at the
last moment via :func:`chrom_for_contigs`.
"""

from __future__ import annotations

import math

VALID_CHROMS: set[str] = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

# Numeric encodings some tools use for sex / mito chromosomes.
_NUMERIC_CHROM_MAP: dict[str, str] = {"23": "X", "24": "Y", "25": "MT", "26": "MT"}

# Common mitochondrial aliases -> canonical "MT".
_MITO_ALIASES: dict[str, str] = {
    "M": "MT",
    "MT": "MT",
    "CHRM": "MT",
    "CHRMT": "MT",
}


def normalise_chrom(raw_chrom: object) -> tuple[str | None, str | None]:
    """Convert a raw chromosome value to a canonical bare chrom string.

    Handles NaN/None/empty, any case of a ``chr`` prefix, bare integers, numeric
    sex/mito encodings (23->X, 24->Y, 25/26->MT), and mito aliases (M->MT).

    Returns ``(bare_chrom, error)``. If ``error`` is not None the value could not
    be normalised and ``bare_chrom`` is None (the variant should be skipped).
    """
    if raw_chrom is None:
        return None, "chromosome is None"
    if isinstance(raw_chrom, float) and math.isnan(raw_chrom):
        return None, "chromosome is NaN (missing value)"

    # Pandas types any chromosome column containing a NaN as float, turning
    # "17" into 17.0 (and numpy.float64 subclasses float). Recover the integer
    # form so numeric chromosomes normalise instead of being rejected as "17.0".
    if isinstance(raw_chrom, float):
        raw_chrom = int(raw_chrom)

    chrom = str(raw_chrom).strip()
    if not chrom or chrom.lower() in ("nan", "none", ".", ""):
        return None, f"chromosome is empty or missing (got {raw_chrom!r})"

    # Same artifact surviving as a string, e.g. "17.0" read from a text cell.
    if chrom.endswith(".0") and chrom[:-2].isdigit():
        chrom = chrom[:-2]

    upper = chrom.upper()

    # Mitochondrial aliases first (before chr-stripping, to catch 'chrM').
    if upper in _MITO_ALIASES:
        return _MITO_ALIASES[upper], None

    # Case-insensitive chr-prefix stripping; preserve case of the remainder.
    bare = chrom[3:] if upper.startswith("CHR") else chrom

    # Numeric alternative encodings (23->X, 24->Y, 25/26->MT).
    bare = _NUMERIC_CHROM_MAP.get(bare, bare)

    # Uppercase for X/Y/MT consistency.
    if bare.upper() in ("X", "Y", "MT"):
        bare = bare.upper()

    if bare not in VALID_CHROMS:
        return None, (
            f"unrecognised chromosome value {raw_chrom!r} "
            f"(normalised to {bare!r}, not in {sorted(VALID_CHROMS)})"
        )
    return bare, None


def contigs_have_chr(contigs) -> bool:
    """Return True if any contig in an iterable is ``chr``-prefixed.

    Probes ``chr1``/``1`` first, then falls back through ``chr2``-``chr5``.
    Defaults to True (assume prefixed) when undeterminable.
    """
    refs = set(contigs)
    if "chr1" in refs:
        return True
    if "1" in refs:
        return False
    for i in range(2, 6):
        if f"chr{i}" in refs:
            return True
        if str(i) in refs:
            return False
    return True


def chrom_for_contigs(bare: str, has_chr: bool) -> str:
    """Convert a bare chromosome to the notation used by a FASTA/VCF."""
    return f"chr{bare}" if has_chr else bare


def detect_series_chr_style(values) -> bool | None:
    """Inspect an iterable of raw chromosome values for reporting.

    Returns True (chr-prefixed), False (bare), or None (unknown/mixed/empty).
    """
    for val in values:
        if val is None:
            continue
        if isinstance(val, float) and math.isnan(val):
            continue
        v = str(val).strip()
        if not v:
            continue
        if v.lower().startswith("chr"):
            return True
        if v in VALID_CHROMS or v in _NUMERIC_CHROM_MAP:
            return False
    return None
