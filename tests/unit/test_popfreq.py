import pytest

from vflank.core.popfreq import (
    GnomadStore,
    example_filename,
    kinds_for,
    parse_common_snp_positions,
    resolve_vcf_for_chrom,
)
from vflank.errors import PopFreqError


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


def test_kinds_for():
    assert kinds_for("genome") == ("genome",)
    assert kinds_for("exome") == ("exome",)
    assert kinds_for("both") == ("genome", "exome")
    with pytest.raises(ValueError):
        kinds_for("nonsense")


def test_resolve_vcf_for_chrom_genome_and_exome(tmp_path):
    # hg19 patterns: gnomad.{genomes,exomes}.r2.1.1.sites.{chrom}.vcf.bgz
    (tmp_path / "gnomad.genomes.r2.1.1.sites.17.vcf.bgz").touch()
    (tmp_path / "gnomad.exomes.r2.1.1.sites.17.vcf.bgz").touch()
    g = resolve_vcf_for_chrom(tmp_path, "17", "hg19", "genome")
    e = resolve_vcf_for_chrom(tmp_path, "17", "hg19", "exome")
    assert g and g.name.startswith("gnomad.genomes")
    assert e and e.name.startswith("gnomad.exomes")
    # missing chrom -> None
    assert resolve_vcf_for_chrom(tmp_path, "5", "hg19", "genome") is None


def test_example_filename():
    assert example_filename("hg19", "exome") == "gnomad.exomes.r2.1.1.sites.1.vcf.bgz"
    assert example_filename("hg38", "genome") == "gnomad.genomes.v4.1.sites.chr1.vcf.bgz"


def test_preflight_errors_when_kind_absent(tmp_path):
    # Only a genome file present; requesting exome must fail fast (not silently
    # fall back to genome-only).
    (tmp_path / "gnomad.genomes.r2.1.1.sites.17.vcf.bgz").touch()
    GnomadStore(tmp_path, "hg19", "genome").preflight(["17"])  # ok
    with pytest.raises(PopFreqError, match="exome"):
        GnomadStore(tmp_path, "hg19", "exome").preflight(["17"])
    with pytest.raises(PopFreqError, match="exome"):
        GnomadStore(tmp_path, "hg19", "both").preflight(["17"])
