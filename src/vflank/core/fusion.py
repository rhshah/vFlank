"""Fusion / SV junction construction from breakpoint pairs.

Builds the chimeric junction sequence for a fusion so a ddPCR probe can span it.

Strand convention (matches iCallSV ``dellyVcf2Tab.py``): ``0`` = plus/reference,
``1`` = minus/complement. The fused product reads 5'->3' as ``partner1 +
partner2`` with **no separator**; partner 1 *ends* at the junction, partner 2
*starts* at it. See docs/research/sv-vcf-input.md for the full derivation
(validated against the unambiguous deletion case).
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(slots=True)
class Breakpoint:
    """A single SV breakpoint. ``pos`` is 1-based; ``strand`` is 0 (+) or 1 (−)."""

    chrom: str
    pos: int
    strand: int


@dataclass(slots=True)
class Fusion:
    """Two breakpoints forming a junction, with optional labels."""

    bp1: Breakpoint
    bp2: Breakpoint
    name: str = ""
    sample: str = ""


@dataclass(slots=True)
class JunctionResult:
    sequence: str
    junction_index: int  # 0-based index where partner 2 begins (= len(partner1))


def _segment(reference, bp: Breakpoint, flank: int, *, donor: bool) -> str:
    """One partner segment of up to ``flank`` bases, oriented to meet the junction.

    ``donor=True``  -> partner 1: its 3' end sits at the junction (ends there).
    ``donor=False`` -> partner 2: its 5' end sits at the junction (starts there).

    Coordinate map (1-based inclusive -> 0-based half-open pysam fetch):
      ``ref[pos-L+1 .. pos]`` = ``fetch(pos-L, pos)``;
      ``ref[pos .. pos+L-1]`` = ``fetch(pos-1, pos-1+L)``.
    """
    pos, L = bp.pos, flank
    if donor:
        if bp.strand == 0:  # plus strand, ending at pos
            seq = reference.fetch(bp.chrom, pos - L, pos)
        else:               # minus strand: revcomp of the bases at/after pos
            seq = reverse_complement(reference.fetch(bp.chrom, pos - 1, pos - 1 + L))
    else:
        if bp.strand == 0:  # revcomp of the bases at/before pos
            seq = reverse_complement(reference.fetch(bp.chrom, pos - L, pos))
        else:               # plus strand, starting at pos
            seq = reference.fetch(bp.chrom, pos - 1, pos - 1 + L)
    return seq.upper()


def build_junction(reference, fusion: Fusion, flank: int) -> JunctionResult:
    """Construct the fusion junction sequence (partner1 + partner2).

    ``flank`` is the bases taken from each partner, so the junction is up to
    ``2*flank`` bp (shorter if a partner runs off a contig end). The probe is
    designed to span ``junction_index``.
    """
    p1 = _segment(reference, fusion.bp1, flank, donor=True)
    p2 = _segment(reference, fusion.bp2, flank, donor=False)
    return JunctionResult(sequence=p1 + p2, junction_index=len(p1))
