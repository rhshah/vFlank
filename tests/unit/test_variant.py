from vflank.core.variant import Variant, validate_allele, validate_coordinates


def test_validate_allele():
    assert validate_allele("ACGT")
    assert validate_allele("-")
    assert validate_allele(".")
    assert validate_allele("")
    assert not validate_allele("ACXT")
    assert not validate_allele("hello")


def test_validate_coordinates():
    assert validate_coordinates(10, 12) is None
    assert validate_coordinates(0, 5) is not None  # start < 1
    assert validate_coordinates(20, 10) is not None  # end < start


def test_variant_raw_chrom_defaults_to_chrom():
    v = Variant(chrom="7", start=1, end=1, ref="A", alt="T")
    assert v.raw_chrom == "7"
    v2 = Variant(chrom="7", start=1, end=1, ref="A", alt="T", raw_chrom="chr7")
    assert v2.raw_chrom == "chr7"
