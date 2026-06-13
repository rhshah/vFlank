"""Flank extraction and masking.

A :class:`FlankSource` is the strategy seam for *where each flank base comes
from*. The reference-backed source implemented here covers modes A (reference)
and B (reference + population mask). Mode C/D (BAM consensus) will add a
``ConsensusFlankSource`` implementing the same protocol.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from .variant import Variant


def mask_sequence(seq: str, region_start_0based: int, positions_1based: list[int]) -> str:
    """Replace bases at the given 1-based genomic positions with 'N'.

    ``seq[0]`` corresponds to genomic position ``region_start_0based + 1``.
    """
    if not positions_1based:
        return seq
    chars = list(seq)
    for pos in positions_1based:
        idx = pos - region_start_0based - 1
        if 0 <= idx < len(chars):
            chars[idx] = "N"
    return "".join(chars)


@dataclass(slots=True)
class FlankResult:
    """Left/right flanks, raw and masked, for one variant."""

    left: str
    right: str
    masked_left: str
    masked_right: str
    covered: int | None = None   # flank positions at/above min_depth (BAM consensus)
    total: int | None = None     # total flank positions (BAM consensus)

    @property
    def n_masked(self) -> int:
        return self.masked_left.count("N") + self.masked_right.count("N")

    @property
    def n_corrected(self) -> int:
        """Flank positions where the masked seq is a real base differing from raw."""
        return sum(
            m != r and m != "N"
            for raw, msk in ((self.left, self.masked_left), (self.right, self.masked_right))
            for r, m in zip(raw, msk, strict=False)
        )

    def upper(self) -> FlankResult:
        """Return an uppercased copy (presentation convenience for callers)."""
        return FlankResult(
            self.left.upper(),
            self.right.upper(),
            self.masked_left.upper(),
            self.masked_right.upper(),
            self.covered,
            self.total,
        )


class FlankSource(Protocol):
    """Strategy for producing flanks around a variant."""

    def fetch(self, variant: Variant) -> FlankResult: ...


class ReferenceFlankSource:
    """Flanks pulled from the reference FASTA, optionally masking common SNPs.

    MAF coordinates are 1-based fully-closed ``[start, end]``; pysam uses 0-based
    half-open ``[start, end)``. The left flank is the ``flank`` bases ending just
    before the variant; the right flank is the ``flank`` bases starting just
    after it. The variant interval itself is excluded from both flanks.
    """

    def __init__(self, reference, gnomad=None, *, flank: int = 200, af_threshold: float = 0.001):
        self.reference = reference
        self.gnomad = gnomad
        self.flank = flank
        self.af_threshold = af_threshold

    def fetch(self, variant: Variant) -> FlankResult:
        left_start_0 = max(0, variant.start - self.flank - 1)
        left_end_0 = variant.start - 1
        right_start_0 = variant.end
        right_end_0 = variant.end + self.flank

        left = self.reference.fetch(variant.chrom, left_start_0, left_end_0)
        right = self.reference.fetch(variant.chrom, right_start_0, right_end_0)

        if self.gnomad is None:
            return FlankResult(left, right, left, right)

        left_snps = self.gnomad.get_positions(
            variant.chrom, left_start_0, left_end_0, self.af_threshold
        )
        right_snps = self.gnomad.get_positions(
            variant.chrom, right_start_0, right_end_0, self.af_threshold
        )
        return FlankResult(
            left,
            right,
            mask_sequence(left, left_start_0, left_snps),
            mask_sequence(right, right_start_0, right_snps),
        )
