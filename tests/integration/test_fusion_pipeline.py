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


def test_fusion_run_with_bam_consensus(tmp_path):
    # A het in a reverse-complemented partner flank must survive as N in the
    # junction (consensus is built in genomic space, then reverse-complemented).
    ref_seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{ref_seq}\n")
    pysam.faidx(str(fasta))

    rs, rl, nr = 20, 50, 24          # reads cover 0-based [20, 70)
    het_g0 = 34                      # genomic 0-based -> het (ref 'G' vs 'T')
    bam = tmp_path / "S1.bam"
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        for i in range(nr):
            seq = list(ref_seq[rs:rs + rl])
            if i >= nr // 2:
                seq[het_g0 - rs] = "T"
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = rs
            a.mapping_quality = 60
            a.cigartuples = [(0, rl)]
            a.query_qualities = pysam.qualitystring_to_array("I" * rl)
            out.write(a)
    pysam.index(str(bam))

    sv = tmp_path / "bp.txt"
    sv.write_text(
        "name\tchr1\tpos1\tstr1\tchr2\tpos2\tstr2\tsample\n"
        "FUSE\t1\t30\t1\t1\t60\t0\tS1\n"   # bp1 str1 (donor revcomp); het at g35 is in its window
    )
    bammap = tmp_path / "map.tsv"
    bammap.write_text(f"S1\t{bam}\n")
    out_fa = tmp_path / "j.fasta"

    result = runner.invoke(app, [
        "fusion", "run", str(sv), "-r", str(fasta), "-g", "hg38",
        "--bam-map", str(bammap), "-f", "10", "--output", str(out_fa),
    ])
    assert result.exit_code == 0, result.output
    lines = out_fa.read_text().splitlines()
    assert lines[0].startswith(">S1__FUSE")          # sample in header
    assert lines[2].startswith(">Masked__S1__FUSE")
    assert lines[1].count("N") == 0                  # raw junction = reference
    assert lines[3].count("N") == 1                  # het -> exactly one N (survived revcomp)
    assert "BAM consensus" in result.output
