"""The :class:`Variant` value object and per-variant validation.

Decoupling the pipeline from pandas rows (the original script iterated
``df.iterrows()`` directly) makes the core logic testable with plain objects.
"""

from __future__ import annotations

from dataclasses import dataclass

VALID_BASES: set[str] = set("ATGCatgc-N")


@dataclass(slots=True)
class Variant:
    """A single small variant in canonical internal form.

    Coordinates are 1-based, fully-closed ``[start, end]`` (MAF convention).
    ``chrom`` is the *bare* chromosome; ``raw_chrom`` preserves the original
    display notation for headers/messages only.
    """

    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    gene: str = "UNKNOWN"
    protein: str = ""
    cdna: str = ""
    sample: str = "SAMPLE"
    raw_chrom: str = ""

    def __post_init__(self) -> None:
        if not self.raw_chrom:
            self.raw_chrom = self.chrom


def validate_allele(allele: str) -> bool:
    """True if an allele is a valid base string (or the '-'/'.' placeholders)."""
    if not allele or allele in (".", "-"):
        return True
    return all(b in VALID_BASES for b in allele)


def validate_coordinates(start: int, end: int) -> str | None:
    """Return an error message if coordinates are invalid, else None."""
    if start < 1:
        return f"start={start} is < 1"
    if end < start:
        return f"end={end} < start={start}"
    return None
