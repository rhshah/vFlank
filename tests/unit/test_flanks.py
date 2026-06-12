from vflank.core.flanks import FlankResult, ReferenceFlankSource, mask_sequence
from vflank.core.variant import Variant


def test_mask_sequence_basic():
    # region starts at 0-based 10, so seq[0] == genomic position 11.
    seq = "ACGTACGTAC"
    # mask genomic position 11 (index 0) and 15 (index 4)
    assert mask_sequence(seq, 10, [11, 15]) == "NCGTNCGTAC"


def test_mask_sequence_out_of_range_ignored():
    assert mask_sequence("ACGT", 10, [5, 100]) == "ACGT"
    assert mask_sequence("ACGT", 10, []) == "ACGT"


class _FakeRef:
    """0-based half-open fetch over an in-memory sequence (single contig)."""

    def __init__(self, seq):
        self.seq = seq

    def fetch(self, chrom, start, end):
        return self.seq[start:end]


def test_reference_flank_source_snp_coordinates():
    # genomic 1..30 = positions; SNP at 1-based pos 15 (index 14)
    seq = "".join("ACGT"[i % 4] for i in range(30))
    ref = _FakeRef(seq)
    v = Variant(chrom="1", start=15, end=15, ref="A", alt="T")
    src = ReferenceFlankSource(ref, gnomad=None, flank=5)
    fr = src.fetch(v)
    # left = bases 10..14 (0-based 9:14), right = bases 16..20 (0-based 15:20)
    assert fr.left == seq[9:14]
    assert fr.right == seq[15:20]
    assert len(fr.left) == 5 and len(fr.right) == 5
    # no gnomad -> masked equals raw
    assert fr.masked_left == fr.left and fr.n_masked == 0


def test_flank_result_n_masked():
    fr = FlankResult("ACGT", "ACGT", "ANGT", "NNGT")
    assert fr.n_masked == 3


def test_flank_result_upper():
    fr = FlankResult("acgt", "tgca", "angt", "nnca").upper()
    assert fr.left == "ACGT" and fr.right == "TGCA"
    assert fr.masked_left == "ANGT" and fr.masked_right == "NNCA"
