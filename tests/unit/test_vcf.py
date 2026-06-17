"""Unit tests for the VCF reader: coordinate mapping, extension detection,
sites-only loading, multi-allelic expansion, symbolic skipping, and best-effort
VEP/SnpEff annotation extraction."""

from __future__ import annotations

import pytest

from vflank.errors import VcfError
from vflank.io.maf import (
    MAF_ALT,
    MAF_CDNA,
    MAF_CHR,
    MAF_END,
    MAF_GENE,
    MAF_PROT,
    MAF_REF,
    MAF_START,
)
from vflank.io.vcf import is_vcf_path, load_vcf, vcf_to_maf_coords

pysam = pytest.importorskip("pysam")


# --- pure coordinate mapping -------------------------------------------------

@pytest.mark.parametrize(
    "pos,ref,alt,expected",
    [
        # SNP: start == end == POS
        (100, "A", "T", (100, 100, "A", "T")),
        # MNV (equal length): a block substitution spanning both bases
        (100, "AC", "GT", (100, 101, "AC", "GT")),
        # insertion (REF prefix of ALT): anchored between POS and POS+1, ref '-'
        (50, "A", "ATG", (50, 51, "-", "TG")),
        # deletion (ALT prefix of REF): deleted span starts after the anchor
        (60, "ATG", "A", (61, 62, "TG", "-")),
        # multi-base anchored insertion
        (10, "AC", "ACGGG", (11, 12, "-", "GGG")),
        # multi-base anchored deletion
        (10, "ACGGG", "AC", (12, 14, "GGG", "-")),
        # complex / non-anchored -> block substitution at POS
        (100, "AT", "GCC", (100, 101, "AT", "GCC")),
        # case-insensitive input is normalised to uppercase
        (5, "a", "g", (5, 5, "A", "G")),
    ],
)
def test_vcf_to_maf_coords(pos, ref, alt, expected):
    assert vcf_to_maf_coords(pos, ref, alt) == expected


def test_is_vcf_path():
    assert is_vcf_path("a.vcf")
    assert is_vcf_path("a.VCF.GZ")
    assert is_vcf_path("/x/y.vcf.bgz")
    assert is_vcf_path("cohort.bcf")
    from pathlib import Path
    assert is_vcf_path(Path("z.vcf.gz"))
    assert not is_vcf_path("variants.maf")
    assert not is_vcf_path("a.tsv")
    # open buffers are never VCF (pysam needs a real path)
    import io
    assert not is_vcf_path(io.StringIO("..."))


# --- load_vcf ----------------------------------------------------------------

_HEADER = "\n".join([
    "##fileformat=VCFv4.2",
    "##contig=<ID=7>",
    '##ALT=<ID=DEL,Description="Deletion">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
])


def _write_vcf(tmp_path, body: str, name: str = "v.vcf"):
    path = tmp_path / name
    path.write_text(_HEADER + "\n" + body + ("\n" if not body.endswith("\n") else ""))
    return path


def test_load_vcf_snp_indels_and_multiallelic(tmp_path):
    vcf = _write_vcf(tmp_path, "\n".join([
        "7\t140453136\t.\tA\tT\t.\t.\t.",      # SNP
        "7\t50\t.\tA\tATG\t.\t.\t.",            # insertion
        "7\t60\t.\tATG\tA\t.\t.\t.",            # deletion
        "7\t100\t.\tA\tG,C\t.\t.\t.",           # multi-allelic -> two rows
    ]))
    df = load_vcf(vcf)

    assert list(df.columns) == [
        MAF_CHR, MAF_START, MAF_END, MAF_REF, MAF_ALT, MAF_GENE, MAF_PROT, MAF_CDNA,
        "Tumor_Sample_Barcode",
    ]
    # 4 records, one multi-allelic with 2 ALTs -> 5 rows.
    assert len(df) == 5
    assert all(df[MAF_CHR] == "7")
    # SNP
    snp = df.iloc[0]
    assert (snp[MAF_START], snp[MAF_END]) == (140453136, 140453136)
    assert (snp[MAF_REF], snp[MAF_ALT]) == ("A", "T")
    # insertion -> ref '-', alt 'TG', anchored [50, 51]
    ins = df.iloc[1]
    assert (ins[MAF_START], ins[MAF_END], ins[MAF_REF], ins[MAF_ALT]) == (50, 51, "-", "TG")
    # deletion -> ref 'TG', alt '-', [61, 62]
    dele = df.iloc[2]
    assert (dele[MAF_START], dele[MAF_END], dele[MAF_REF], dele[MAF_ALT]) == (61, 62, "TG", "-")
    # multi-allelic A>G and A>C both at 100
    assert list(df.iloc[3:5][MAF_ALT]) == ["G", "C"]
    assert all(df.iloc[3:5][MAF_START] == 100)
    # sites-only -> every row gets the placeholder sample
    assert all(df["Tumor_Sample_Barcode"] == "SAMPLE")
    # no annotation -> gene UNKNOWN, HGVS blank
    assert all(df[MAF_GENE] == "UNKNOWN")
    assert all(df[MAF_PROT] == "")


def test_load_vcf_skips_symbolic_and_errors_when_empty(tmp_path):
    # Only a symbolic SV ALT -> nothing usable -> typed error.
    vcf = _write_vcf(tmp_path, "7\t200\t.\tA\t<DEL>\t.\t.\tEND=300")
    with pytest.raises(VcfError, match="No small variants"):
        load_vcf(vcf)


def test_load_vcf_open_error(tmp_path):
    missing = tmp_path / "nope.vcf"
    with pytest.raises(VcfError, match="Could not open VCF"):
        load_vcf(missing)


def test_load_vcf_n_rows_caps(tmp_path):
    vcf = _write_vcf(tmp_path, "\n".join(
        f"7\t{100 + i}\t.\tA\tT\t.\t.\t." for i in range(5)
    ))
    assert len(load_vcf(vcf, n_rows=2)) == 2


def test_load_vcf_vep_csq_annotation(tmp_path):
    header = "\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=7>",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
        'from Ensembl VEP. Format: Allele|Consequence|SYMBOL|HGVSc|HGVSp">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ])
    path = tmp_path / "vep.vcf"
    path.write_text(
        header + "\n"
        + "7\t140453136\t.\tA\tT\t.\t.\tCSQ=T|missense_variant|BRAF|c.1799TA|p.Val600Glu\n"
    )
    df = load_vcf(path)
    assert df.iloc[0][MAF_GENE] == "BRAF"
    assert df.iloc[0][MAF_CDNA] == "c.1799TA"
    assert df.iloc[0][MAF_PROT] == "p.Val600Glu"


def test_load_vcf_snpeff_ann_annotation(tmp_path):
    header = "\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=7>",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Gene_Name | HGVS.c | HGVS.p'\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ])
    path = tmp_path / "snpeff.vcf"
    path.write_text(
        header + "\n"
        + "7\t140453136\t.\tA\tT\t.\t.\tANN=T|missense_variant|BRAF|c.1799TA|p.Val600Glu\n"
    )
    df = load_vcf(path)
    assert df.iloc[0][MAF_GENE] == "BRAF"
    assert df.iloc[0][MAF_CDNA] == "c.1799TA"
    assert df.iloc[0][MAF_PROT] == "p.Val600Glu"
