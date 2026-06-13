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

from .flanks import mask_sequence

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
    masked_sequence: str
    junction_index: int  # 0-based index where partner 2 begins (= len(partner1))
    n_masked: int = 0


def _segment(
    reference, bp: Breakpoint, flank: int, *, donor: bool, gnomad=None,
    af_threshold: float = 0.001, bam_source=None,
) -> tuple[str, str]:
    """Return ``(raw, masked)`` for one partner segment, oriented to the junction.

    ``donor=True``  -> partner 1: its 3' end sits at the junction (ends there).
    ``donor=False`` -> partner 2: its 5' end sits at the junction (starts there).

    Strand 0 takes ``ref[pos-L+1 .. pos]``; strand 1 takes ``ref[pos .. pos+L-1]``;
    whether the partner is reverse-complemented depends on its role (see table in
    docs/research/sv-vcf-input.md). Masking is applied in **genomic (plus-strand)
    space before any reverse-complement** — since ``revcomp(N) == N`` this equals
    masking the oriented segment, and reuses ``get_positions`` + ``mask_sequence``.
    """
    pos, L = bp.pos, flank
    if bp.strand == 0:
        start0, end0 = max(0, pos - L), pos
        rc = not donor
    else:
        start0, end0 = pos - 1, pos - 1 + L
        rc = donor

    raw = reference.fetch(bp.chrom, start0, end0).upper()
    masked = raw
    if bam_source is not None:
        # Patient consensus in genomic (plus-strand) space; revcomp below.
        gnomad_pos = (
            set(gnomad.get_positions(bp.chrom, start0, end0, af_threshold))
            if gnomad is not None else set()
        )
        masked, _covered = bam_source.consensus(bp.chrom, start0, end0, raw, gnomad_pos)
    elif gnomad is not None:
        snps = gnomad.get_positions(bp.chrom, start0, end0, af_threshold)
        masked = mask_sequence(raw, start0, snps)
    if rc:
        raw, masked = reverse_complement(raw), reverse_complement(masked)
    return raw, masked


def build_junction(
    reference, fusion: Fusion, flank: int, gnomad=None, af_threshold: float = 0.001,
    bam_source=None,
) -> JunctionResult:
    """Construct the fusion junction sequence (partner1 + partner2).

    ``flank`` is the bases taken from each partner, so the junction is up to
    ``2*flank`` bp (shorter if a partner runs off a contig end). The probe is
    designed to span ``junction_index``. ``masked_sequence`` is the gnomAD-masked
    junction, or — when ``bam_source`` is given — the per-sample patient consensus
    (built in genomic space before reverse-complement).
    """
    raw1, masked1 = _segment(reference, fusion.bp1, flank, donor=True,
                             gnomad=gnomad, af_threshold=af_threshold, bam_source=bam_source)
    raw2, masked2 = _segment(reference, fusion.bp2, flank, donor=False,
                             gnomad=gnomad, af_threshold=af_threshold, bam_source=bam_source)
    sequence = raw1 + raw2
    masked_sequence = masked1 + masked2
    return JunctionResult(
        sequence=sequence,
        masked_sequence=masked_sequence,
        junction_index=len(raw1),
        n_masked=masked_sequence.count("N") - sequence.count("N"),
    )
