from vflank.core.consensus import (
    ConsensusPolicy,
    _parse_consensus_fasta,
    apply_lowcov_overlay,
)


def test_parse_consensus_fasta():
    assert _parse_consensus_fasta(">chr1\nACGT\nACG\n") == "ACGTACG"
    assert _parse_consensus_fasta(">x\nacgt\n") == "ACGT"


def test_overlay_keeps_adequate_depth():
    # Every position has depth >= min_depth -> the samtools call is kept, even N.
    p = ConsensusPolicy(min_depth=20)
    assert apply_lowcov_overlay("ACNT", [30, 30, 30, 30], 100, "ACGT", set(), p) == "ACNT"


def test_overlay_lowcov_n():
    p = ConsensusPolicy(min_depth=20, lowcov="n")
    # index 1 is low-coverage -> N regardless of reference
    assert apply_lowcov_overlay("ACGT", [30, 5, 30, 30], 100, "ATGT", set(), p) == "ANGT"


def test_overlay_lowcov_reference():
    p = ConsensusPolicy(min_depth=20, lowcov="reference")
    # index 1 low -> reference base ('T'); others kept
    assert apply_lowcov_overlay("ACGT", [30, 5, 30, 30], 100, "ATGT", set(), p) == "ATGT"


def test_overlay_lowcov_gnomad():
    p = ConsensusPolicy(min_depth=20, lowcov="gnomad")
    # region_start0=100 -> idx0=pos101 (in gnomad -> N), idx1=pos102 (not -> ref 'G')
    out = apply_lowcov_overlay("ACGT", [5, 5, 30, 30], 100, "GGGT", {101}, p)
    assert out == "NGGT"
