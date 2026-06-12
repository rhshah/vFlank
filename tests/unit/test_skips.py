from vflank.core.skips import categorize_skip


def test_categories():
    chrom_cat = "chromosome missing/invalid"
    assert categorize_skip("TP53 — chromosome is NaN (missing value)") == chrom_cat
    assert categorize_skip("unrecognised chromosome value '99'") == chrom_cat
    assert categorize_skip("GENE: non-numeric position (start='x')") == "non-numeric position"
    assert categorize_skip("GENE 7:10 — invalid allele (ref='Q')") == "invalid allele"
    assert categorize_skip("GENE 7:10 — fetch error: boom") == "fetch error"
    assert categorize_skip("GENE 7:10 — end=5 < start=10") == "invalid coordinates"
    assert categorize_skip("GENE 7:0 — start=0 is < 1") == "invalid coordinates"
    assert categorize_skip("something unexpected") == "other"
