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
