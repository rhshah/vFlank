"""Synthetic-BAM tests for the samtools-consensus engine (Phase 1)."""

from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from vflank.core.consensus import (  # noqa: E402
    ConsensusPolicy,
    run_samtools_consensus,
    window_depth,
)

BASE = "ACGTG"
READ_LEN = 30
REF_START = 10  # 0-based; reads cover 0-based [10, 40)


def _template():
    return [BASE[i % len(BASE)] for i in range(READ_LEN)]


def _make_bam(tmp_path, *, n_reads=20, het_idx=None, het_alt="T"):
    """A BAM of identical reads, optionally heterozygous (50/50) at one index."""
    bam = tmp_path / "t.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 60}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(n_reads):
            seq = _template()
            if het_idx is not None and i >= n_reads // 2:
                seq[het_idx] = het_alt if seq[het_idx] != het_alt else "A"
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = REF_START
            a.mapping_quality = 60
            a.cigartuples = [(0, READ_LEN)]
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            out.write(a)
    pysam.index(str(bam))
    return str(bam)


def _genomic_to_read_idx(pos_1based):
    return pos_1based - REF_START - 1  # template index for a genomic position


def test_consensus_homozygous_matches_reads(tmp_path):
    bam = _make_bam(tmp_path)
    seq = run_samtools_consensus(bam, "chr1", 22, 32, ConsensusPolicy())
    assert len(seq) == 11
    # genomic pos 30 -> template base
    assert seq[30 - 22] == BASE[_genomic_to_read_idx(30) % len(BASE)]


def test_heterozygous_becomes_N(tmp_path):
    # het at genomic 25 (read idx 14: ref 'G' vs alt 'T') -> default consensus N
    bam = _make_bam(tmp_path, het_idx=14, het_alt="T")
    seq = run_samtools_consensus(bam, "chr1", 22, 32, ConsensusPolicy())
    assert seq[25 - 22] == "N"


def test_heterozygous_iupac(tmp_path):
    bam = _make_bam(tmp_path, het_idx=14, het_alt="T")
    seq = run_samtools_consensus(bam, "chr1", 22, 32, ConsensusPolicy(het_char="iupac"))
    assert seq[25 - 22] in "RYSWKM"  # an ambiguity code, not N


def test_window_depth(tmp_path):
    bam = _make_bam(tmp_path, n_reads=20)
    # 0-based [21, 32): all covered by 20 reads
    assert window_depth(bam, "chr1", 21, 32, ConsensusPolicy()) == [20] * 11
    # a region outside the read span -> zero depth
    assert window_depth(bam, "chr1", 45, 50, ConsensusPolicy()) == [0] * 5


def test_consensus_flank_source(tmp_path):
    from vflank.core.consensus import BamConsensusSource, ConsensusFlankSource
    from vflank.core.variant import Variant
    from vflank.io.reference import ReferenceFasta

    ref_seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{ref_seq}\n")
    pysam.faidx(str(fasta))

    rs, rl, nr = 25, 40, 24          # 24 reads cover 1-based 26..65
    het_g0 = 44                      # genomic 0-based 44 -> het (ref 'A' vs 'G')
    homalt_g0 = 54                   # genomic 0-based 54 -> hom-ALT (ref 'G' -> 'T')
    bam = tmp_path / "s.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(nr):
            seq = list(ref_seq[rs:rs + rl])
            seq[homalt_g0 - rs] = "T"             # all reads: hom-ALT
            if i >= nr // 2:
                seq[het_g0 - rs] = "G"            # half reads: het alt
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = rs
            a.mapping_quality = 60
            a.cigartuples = [(0, rl)]
            a.query_qualities = pysam.qualitystring_to_array("I" * rl)
            out.write(a)
    pysam.index(str(bam))

    ref = ReferenceFasta(str(fasta))
    src = ConsensusFlankSource(ref, BamConsensusSource(str(bam)), gnomad=None, flank=10)
    fr = src.fetch(Variant(chrom="1", start=50, end=50, ref="A", alt="T"))
    # left window 1-based 40..49; het at genomic 45 -> idx 5 -> N
    assert fr.masked_left[5] == "N"
    # right window 1-based 51..60; hom-ALT at genomic 55 -> idx 4 -> patient 'T' (ref 'G')
    assert fr.masked_right[4] == "T" and fr.right[4] == "G"
    # everything else equals the reference flank (raw)
    assert fr.masked_left[:5] + fr.masked_left[6:] == fr.left[:5] + fr.left[6:]
    ref.close()
