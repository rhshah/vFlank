from vflank.core.popfreq import parse_common_snp_positions


def _vcf_line(pos, ref, alt, info):
    # CHROM POS ID REF ALT QUAL FILTER INFO
    return f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}"


def test_common_snp_above_threshold():
    rows = [_vcf_line(100, "A", "G", "AF=0.25")]
    assert parse_common_snp_positions(rows, 0.001) == [100]


def test_below_threshold_excluded():
    rows = [_vcf_line(100, "A", "G", "AF=0.0001")]
    assert parse_common_snp_positions(rows, 0.001) == []


def test_indels_excluded():
    rows = [
        _vcf_line(100, "AT", "A", "AF=0.5"),   # deletion: multi-base REF
        _vcf_line(200, "A", "ATG", "AF=0.5"),  # insertion: multi-base ALT
    ]
    assert parse_common_snp_positions(rows, 0.001) == []


def test_multiallelic_uses_max_af():
    rows = [_vcf_line(100, "A", "G,C", "AF=0.0001,0.4")]
    assert parse_common_snp_positions(rows, 0.001) == [100]


def test_grpmax_field_and_nan_guard():
    rows = [
        _vcf_line(100, "A", "G", "AC=1;AF_grpmax=0.2"),
        _vcf_line(200, "A", "G", "AF=nan"),  # NaN must not crash or pass
    ]
    assert parse_common_snp_positions(rows, 0.001) == [100]


def test_star_allele_skipped():
    rows = [_vcf_line(100, "A", "*", "AF=0.9")]
    assert parse_common_snp_positions(rows, 0.001) == []
