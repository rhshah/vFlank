from vflank.core.chrom import (
    chrom_for_contigs,
    contigs_have_chr,
    detect_series_chr_style,
    normalise_chrom,
)


def test_normalise_plain_and_prefixed():
    assert normalise_chrom("7") == ("7", None)
    assert normalise_chrom("chr7") == ("7", None)
    assert normalise_chrom("CHR7") == ("7", None)
    assert normalise_chrom("Chr7") == ("7", None)


def test_normalise_sex_and_mito():
    assert normalise_chrom("X") == ("X", None)
    assert normalise_chrom("chrx") == ("X", None)
    assert normalise_chrom("23") == ("X", None)
    assert normalise_chrom("24") == ("Y", None)
    assert normalise_chrom("25")[0] == "MT"
    assert normalise_chrom("M") == ("MT", None)
    assert normalise_chrom("chrM") == ("MT", None)
    assert normalise_chrom("MT") == ("MT", None)


def test_normalise_float_chromosome_from_pandas():
    # A NaN-containing column makes pandas store "17" as float 17.0.
    assert normalise_chrom(17.0) == ("17", None)
    assert normalise_chrom("17.0") == ("17", None)
    assert normalise_chrom(23.0) == ("X", None)  # float numeric sex encoding
    assert normalise_chrom(24.0) == ("Y", None)


def test_normalise_invalid():
    assert normalise_chrom(None)[0] is None
    assert normalise_chrom(float("nan"))[0] is None
    assert normalise_chrom("")[0] is None
    assert normalise_chrom("banana")[0] is None
    assert normalise_chrom("99")[0] is None


def test_contig_style_detection():
    assert contigs_have_chr(["chr1", "chr2"]) is True
    assert contigs_have_chr(["1", "2"]) is False
    assert contigs_have_chr(["chr3", "chr4"]) is True  # falls through past chr1
    assert chrom_for_contigs("7", True) == "chr7"
    assert chrom_for_contigs("7", False) == "7"


def test_series_style():
    assert detect_series_chr_style(["chr1", "chr2"]) is True
    assert detect_series_chr_style(["1", "2"]) is False
    assert detect_series_chr_style([None, "", float("nan")]) is None
