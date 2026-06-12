"""End-to-end test: build a tiny indexed FASTA and run `vflank fusion run`."""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

pysam = pytest.importorskip("pysam")

from vflank.cli.app import app  # noqa: E402
from vflank.core.fusion import reverse_complement  # noqa: E402

runner = CliRunner()


def test_fusion_run_deletion_like(tmp_path):
    # Two contigs so we can check an inter-chromosomal junction.
    seq1 = "".join("ACGTG"[i % 5] for i in range(60))
    seq2 = "".join("TTAGC"[i % 5] for i in range(60))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{seq1}\n>2\n{seq2}\n")
    pysam.faidx(str(fasta))

    # Columns deliberately scrambled + extra column, to exercise name-driven parsing.
    sv = tmp_path / "bp.txt"
    sv.write_text(
        "name\tstr2\tchr1\tpos1\tstr1\tchr2\tpos2\textra\n"
        "FUSE\t1\t1\t20\t0\t2\t30\tjunk\n"
    )
    out = tmp_path / "j.fasta"

    result = runner.invoke(app, [
        "fusion", "run", str(sv), "--ref-genome", str(fasta),
        "--genome-build", "hg38", "--flank", "5", "--output", str(out),
    ])
    assert result.exit_code == 0, result.output

    lines = out.read_text().splitlines()
    assert lines[0].startswith(">FUSE__1_20_0__2_30_1__j5")
    # partner1 (chr1, str0): plus ending at pos1 -> seq1[15:20]
    # partner2 (chr2, str1): plus starting at pos2 -> seq2[29:34]
    assert lines[1] == (seq1[15:20] + seq2[29:34]).upper()


def test_fusion_run_minus_strand_revcomp(tmp_path):
    seq1 = "".join("ACGTG"[i % 5] for i in range(60))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{seq1}\n")
    pysam.faidx(str(fasta))

    sv = tmp_path / "bp.txt"
    sv.write_text(
        "chr1\tpos1\tstr1\tchr2\tpos2\tstr2\n"
        "1\t25\t1\t1\t40\t0\n"  # both revcomp branches
    )
    out = tmp_path / "j.fasta"
    result = runner.invoke(app, [
        "fusion", "run", str(sv), "-r", str(fasta), "-g", "hg38", "-f", "5", "-o", str(out),
    ])
    assert result.exit_code == 0, result.output
    seq = out.read_text().splitlines()[1]
    expected = reverse_complement(seq1[24:29]).upper() + reverse_complement(seq1[35:40]).upper()
    assert seq == expected
