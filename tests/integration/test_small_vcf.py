"""End-to-end: `vflank small run` auto-detects a VCF input and produces the
same FASTA shape as the MAF path (variant shown literally, flanks from the
reference)."""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

pysam = pytest.importorskip("pysam")

from vflank.cli.app import app  # noqa: E402

runner = CliRunner()

_VCF_HEADER = "\n".join([
    "##fileformat=VCFv4.2",
    "##contig=<ID=7>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
])


def _write_reference(tmp_path):
    # 80 bp single contig named '7' (bare notation, like the MAF tests).
    seq = "".join("ACGT"[i % 4] for i in range(80))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">7\n{seq}\n")
    pysam.faidx(str(fasta))
    return fasta, seq


def _write_vcf(tmp_path, body: str, name: str = "v.vcf"):
    path = tmp_path / name
    path.write_text(_VCF_HEADER + "\n" + body + "\n")
    return path


def test_run_vcf_snp_matches_reference_flanks(tmp_path):
    fasta, seq = _write_reference(tmp_path)
    vcf = _write_vcf(tmp_path, "7\t40\t.\tA\tG\t.\t.\t.")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(vcf), "--ref-genome", str(fasta),
        "--genome-build", "hg19", "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    # CLI announces the VCF path was taken.
    assert "Loading VCF" in result.output

    lines = out.read_text().splitlines()
    assert len(lines) == 4
    # No annotation -> gene UNKNOWN, blank HGVS; key is CHROM_POS_REF_ALT.
    assert lines[0] == ">UNKNOWN__.__.__7_40_A_G"
    assert lines[1] == f"{seq[34:39].upper()}[A/G]{seq[40:45].upper()}"
    assert lines[2] == ">Masked__UNKNOWN__.__.__7_40_A_G"


def test_run_vcf_deletion_coordinates(tmp_path):
    fasta, seq = _write_reference(tmp_path)
    # ATG->A deletion anchored at POS 40 -> deleted span [41, 42], ref 'TG', alt '-'.
    vcf = _write_vcf(tmp_path, "7\t40\t.\tATG\tA\t.\t.\t.")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(vcf), "--ref-genome", str(fasta),
        "-g", "hg19", "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    lines = out.read_text().splitlines()
    assert lines[0] == ">UNKNOWN__.__.__7_41_TG_-"
    # left = ref[35:40] (0-based, ends before start-1=40), right = ref[42:47].
    assert lines[1] == f"{seq[35:40].upper()}[TG/-]{seq[42:47].upper()}"


def test_run_vcf_insertion_coordinates(tmp_path):
    fasta, seq = _write_reference(tmp_path)
    # A->ATG insertion anchored at POS 40 -> Start=40, End=41, ref '-', alt 'TG'.
    vcf = _write_vcf(tmp_path, "7\t40\t.\tA\tATG\t.\t.\t.")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(vcf), "--ref-genome", str(fasta),
        "-g", "hg19", "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    lines = out.read_text().splitlines()
    assert lines[0] == ">UNKNOWN__.__.__7_40_-_TG"
    # left = ref[34:39] (ends before start-1=39), right = ref[41:46] (from end=41).
    assert lines[1] == f"{seq[34:39].upper()}[-/TG]{seq[41:46].upper()}"


def test_run_vcf_multiallelic_emits_both(tmp_path):
    fasta, _ = _write_reference(tmp_path)
    vcf = _write_vcf(tmp_path, "7\t40\t.\tA\tG,C\t.\t.\t.")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(vcf), "--ref-genome", str(fasta),
        "-g", "hg19", "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    # Two distinct ALTs -> two variants -> 8 lines, both kept under dedup.
    text = out.read_text()
    assert text.count(">") == 4  # 2 variants x (raw + masked)
    assert "7_40_A_G" in text and "7_40_A_C" in text


def test_run_vcf_samples_filter_ignored(tmp_path):
    fasta, _ = _write_reference(tmp_path)
    vcf = _write_vcf(tmp_path, "7\t40\t.\tA\tG\t.\t.\t.")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(vcf), "--ref-genome", str(fasta),
        "-g", "hg19", "--flank", "5", "--output", str(out), "--samples", "S1",
    ])
    assert result.exit_code == 0, result.output
    assert "ignored for VCF" in result.output
    assert len(out.read_text().splitlines()) == 4  # still processed


def test_inspect_vcf(tmp_path):
    vcf = _write_vcf(tmp_path, "\n".join([
        "7\t40\t.\tA\tG\t.\t.\t.",
        "7\t50\t.\tA\tATG\t.\t.\t.",
    ]))
    result = runner.invoke(app, ["small", "inspect", str(vcf)])
    assert result.exit_code == 0, result.output
    assert "VCF Inspect" in result.output
