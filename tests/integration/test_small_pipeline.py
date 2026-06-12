"""End-to-end test: build a tiny indexed FASTA and run `vflank small run`."""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

pysam = pytest.importorskip("pysam")

from vflank.cli.app import app  # noqa: E402

runner = CliRunner()


def _write_reference(tmp_path):
    # 60 bp single contig named 'chr1' (FASTA uses chr-prefixed notation).
    seq = "".join("ACGT"[i % 4] for i in range(60))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(fasta))  # creates ref.fasta.fai
    return fasta, seq


def _write_maf(tmp_path):
    # MAF uses bare '1' notation to exercise normalisation against a chr FASTA.
    header = "\t".join([
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ])
    row = "\t".join(["TP53", "1", "30", "30", "A", "T", "SAMPLE_1"])
    maf = tmp_path / "variants.maf"
    maf.write_text(header + "\n" + row + "\n")
    return maf


def test_run_produces_expected_fasta(tmp_path):
    fasta, seq = _write_reference(tmp_path)
    maf = _write_maf(tmp_path)
    out = tmp_path / "out.fasta"

    report = tmp_path / "report.tsv"
    result = runner.invoke(app, [
        "small", "run", str(maf),
        "--ref-genome", str(fasta),
        "--genome-build", "hg38",
        "--flank", "5",
        "--output", str(out),
        "--report", str(report),
    ])

    assert result.exit_code == 0, result.output
    assert out.exists()

    # The run report is written and records the single processed variant.
    assert report.exists()
    rtext = report.read_text()
    assert "# processed\t1" in rtext
    assert "TP53" in rtext.splitlines()[-1]

    lines = out.read_text().splitlines()
    # Two records (raw + masked) = 4 lines.
    assert len(lines) == 4
    assert lines[0].startswith(">SAMPLE_1__TP53")
    assert lines[2].startswith(">Masked__SAMPLE_1__TP53")

    # Variant at 1-based pos 30 (index 29). flank=5.
    expected_left = seq[24:29]   # positions 25..29
    expected_right = seq[30:35]  # positions 31..35
    assert lines[1] == f"{expected_left}[A/T]{expected_right}"
    # No gnomAD dir -> masked record identical to raw.
    assert lines[3].endswith(f"{expected_left}[A/T]{expected_right}")


def test_truncated_flank_is_emitted_not_dropped(tmp_path):
    # 60 bp contig; variant near the end with a flank that runs off the contig.
    fasta, seq = _write_reference(tmp_path)
    header = "\t".join([
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ])
    # pos 58 of a 60 bp contig: right flank can only be 2 bp, not the requested 10.
    maf = tmp_path / "edge.maf"
    maf.write_text(header + "\n" + "\t".join(["TP53", "1", "58", "58", "A", "T", "S1"]) + "\n")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(maf), "--ref-genome", str(fasta),
        "--flank", "10", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output

    # The record must still be emitted (truncation warns, never silently drops).
    raw = out.read_text().splitlines()[1]
    right = raw.split("]", 1)[1]
    assert len(right) == 2  # only 2 bp available past position 58
    assert raw.split("[", 1)[0] == seq[47:57].upper()  # full 10 bp left flank


def test_bare_named_single_contig_fasta_resolves(tmp_path):
    # A FASTA whose only contig is bare '7' defeats best-effort chr detection
    # (which probes chr1-5). ReferenceFasta.contig must fall back to the form
    # that actually exists rather than KeyError-ing on every variant.
    seq = "".join("ACGT"[i % 4] for i in range(80))
    fasta = tmp_path / "ref7.fasta"
    fasta.write_text(f">7\n{seq}\n")
    pysam.faidx(str(fasta))

    header = "\t".join([
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ])
    maf = tmp_path / "v.maf"
    maf.write_text(header + "\n" + "\t".join(["BRAF", "7", "40", "40", "A", "T", "S1"]) + "\n")
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(maf), "--ref-genome", str(fasta),
        "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    raw = out.read_text().splitlines()[1]
    assert raw == f"{seq[34:39].upper()}[A/T]{seq[40:45].upper()}"


def test_float_typed_chromosome_column_is_processed(tmp_path):
    # When a MAF's Chromosome column contains any blank, pandas types the whole
    # column as float, turning "17" into 17.0. The valid row must still process
    # (regression for a real MSK MAF where this lost every variant).
    seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref17.fasta"
    fasta.write_text(f">17\n{seq}\n")
    pysam.faidx(str(fasta))

    header = "\t".join([
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ])
    maf = tmp_path / "mixed.maf"
    maf.write_text(
        header + "\n"
        + "\t".join(["TP53", "17", "40", "40", "A", "T", "S1"]) + "\n"  # valid
        + "\t".join(["", "", "50", "50", "", "", ""]) + "\n"            # blanks -> float column
    )
    out = tmp_path / "out.fasta"
    result = runner.invoke(app, [
        "small", "run", str(maf), "--ref-genome", str(fasta),
        "--genome-build", "hg19", "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    lines = out.read_text().splitlines()
    assert len(lines) == 4  # exactly one variant processed (raw + masked records)
    assert lines[1] == f"{seq[34:39].upper()}[A/T]{seq[40:45].upper()}"


def test_pop_data_exome_without_files_errors(tmp_path):
    # --pop-data exome but the dir has no exome VCFs -> fail fast, do not
    # silently mask with genome-only (or nothing).
    fasta, _ = _write_reference(tmp_path)  # contig 'chr1'
    header = "\t".join([
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ])
    maf = tmp_path / "v.maf"
    maf.write_text(header + "\n" + "\t".join(["G", "1", "30", "30", "A", "T", "S1"]) + "\n")
    empty_dir = tmp_path / "gnomad"
    empty_dir.mkdir()
    out = tmp_path / "out.fasta"

    result = runner.invoke(app, [
        "small", "run", str(maf), "--ref-genome", str(fasta), "--genome-build", "hg38",
        "--pop-vcf-dir", str(empty_dir), "--pop-data", "exome",
        "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 1
    assert "exome" in result.output


def test_missing_required_column_errors(tmp_path):
    fasta, _ = _write_reference(tmp_path)
    bad = tmp_path / "bad.maf"
    bad.write_text("Hugo_Symbol\tChromosome\nTP53\t1\n")
    out = tmp_path / "out.fasta"
    result = runner.invoke(app, [
        "small", "run", str(bad), "--ref-genome", str(fasta), "--output", str(out),
    ])
    assert result.exit_code == 1
