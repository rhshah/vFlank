from vflank.core.fusion import (
    Breakpoint,
    Fusion,
    build_junction,
    reverse_complement,
)


class _FakeRef:
    """0-based half-open fetch over in-memory per-chromosome sequences."""

    def __init__(self, seqs):
        self.seqs = seqs

    def fetch(self, chrom, start, end):
        return self.seqs[chrom][start:end]


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAC") == "GTTT"
    assert reverse_complement("acgtN") == "Nacgt"


# chrom "1" length 50, deterministic but non-palindromic.
SEQ = "".join("ACGTG"[i % 5] for i in range(50))


def _ref():
    return _FakeRef({"1": SEQ, "2": SEQ[::-1]})


def test_deletion_junction_is_plus_plus():
    # CT=3to5 -> (str1=0, str2=1): truth = plus[..pos1] + plus[pos2..], no revcomp.
    ref = _ref()
    fus = Fusion(Breakpoint("1", 20, 0), Breakpoint("1", 30, 1))
    res = build_junction(ref, fus, flank=5)
    assert res.sequence == (SEQ[15:20] + SEQ[29:34]).upper()
    assert res.junction_index == 5


def test_all_four_strand_combinations():
    ref = _ref()
    L = 5

    def seg_donor(pos, strand):
        return (SEQ[pos-L:pos] if strand == 0
                else reverse_complement(SEQ[pos-1:pos-1+L])).upper()

    def seg_acceptor(pos, strand):
        return (reverse_complement(SEQ[pos-L:pos]) if strand == 0
                else SEQ[pos-1:pos-1+L]).upper()

    for s1 in (0, 1):
        for s2 in (0, 1):
            fus = Fusion(Breakpoint("1", 25, s1), Breakpoint("1", 40, s2))
            res = build_junction(ref, fus, flank=L)
            expected = seg_donor(25, s1) + seg_acceptor(40, s2)
            assert res.sequence == expected, (s1, s2)


def test_interchromosomal_uses_both_contigs():
    ref = _ref()
    fus = Fusion(Breakpoint("1", 20, 0), Breakpoint("2", 30, 1))
    res = build_junction(ref, fus, flank=5)
    assert res.sequence[:5] == SEQ[15:20].upper()           # partner1 from chrom 1
    assert res.sequence[5:] == ref.seqs["2"][29:34].upper()  # partner2 from chrom 2


class _FakeGnomad:
    """get_positions returning 1-based positions that fall in [start0, end0)."""

    def __init__(self, positions):
        self.positions = positions

    def get_positions(self, chrom, start0, end0, af):
        return [p for p in self.positions if start0 < p <= end0]


def test_build_junction_masks_partner_flanks():
    ref = _ref()
    # bp1 donor, str0 (not revcomp): window 0-based [15,20); mask 1-based pos 17.
    fus = Fusion(Breakpoint("1", 20, 0), Breakpoint("1", 40, 1))
    res = build_junction(ref, fus, flank=5, gnomad=_FakeGnomad([17]))
    assert res.masked_sequence[1] == "N"          # idx 17-15-1 = 1
    assert res.sequence[1] != "N"                 # raw is unchanged
    assert res.n_masked == 1


def test_masking_survives_reverse_complement():
    # bp1 donor str1 -> the partner IS reverse-complemented; masking (done in
    # genomic space before revcomp) must still land an N, since revcomp(N)=N.
    ref = _ref()
    fus = Fusion(Breakpoint("1", 20, 1), Breakpoint("1", 40, 1))
    res = build_junction(ref, fus, flank=5, gnomad=_FakeGnomad([22]))
    assert res.masked_sequence.count("N") == 1
    assert res.n_masked == 1
