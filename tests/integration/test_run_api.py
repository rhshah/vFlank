"""End-to-end tests for the library entrypoints run_small / run_fusion.

These are the functions a web service / notebook calls. They also exercise
buffer input (an in-memory StringIO instead of a path), the way an uploaded file
arrives in a service.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from vflank.errors import VflankError  # noqa: E402
from vflank.pipeline import run_fusion, run_small  # noqa: E402


def _reference(tmp_path):
    seq = "".join("ACGT"[i % 4] for i in range(60))
    fa = tmp_path / "ref.fasta"
    fa.write_text(f">1\n{seq}\n")
    pysam.faidx(str(fa))
    return fa, seq


_MAF = (
    "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
    "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n"
    "TP53\t1\t30\t30\tA\tT\tS1\n"
)
_SV = "name\tchr1\tpos1\tstr1\tchr2\tpos2\tstr2\nFUSE\t1\t20\t0\t1\t40\t0\n"


def test_run_small_from_buffer(tmp_path):
    fa, seq = _reference(tmp_path)
    result = run_small(
        io.StringIO(_MAF), genome_build="hg19", ref_source="file", ref_genome=fa, flank=5
    )
    assert result.n_processed == 1
    assert result.n_skipped == 0
    assert len(result.records) == 2                      # raw + masked
    assert result.records[0].startswith(">TP53__.__.__1_30_A_T\n")
    assert result.records[1].startswith(">Masked__TP53__.__.__1_30_A_T\n")
    # No gnomAD here, so masked == raw apart from the header.
    assert result.records[0].splitlines()[1] == result.records[1].splitlines()[1]
    assert result.rows[0]["Gene"] == "TP53"
    assert result.ref_api_requests is None               # file backend, not API


def test_run_small_emit_primer3(tmp_path):
    fa, _ = _reference(tmp_path)
    result = run_small(
        io.StringIO(_MAF), genome_build="hg19", ref_source="file", ref_genome=fa,
        flank=5, emit_primer3=True,
    )
    assert len(result.primer3) == 1
    assert result.primer3[0].id.startswith("TP53__")


def test_run_fusion_from_buffer(tmp_path):
    fa, _ = _reference(tmp_path)
    result = run_fusion(
        io.StringIO(_SV), genome_build="hg19", ref_source="file", ref_genome=fa, flank=5
    )
    assert result.n_processed == 1
    assert len(result.records) == 2
    assert result.records[0].startswith(">FUSE__")


def test_run_small_validates_options_before_io():
    # Bad build -> VflankError before any file/network access (ref_genome unused).
    with pytest.raises(VflankError, match="genome-build"):
        run_small(io.StringIO(_MAF), genome_build="hg99", ref_genome=Path("/does/not/exist"))


def test_run_small_file_backend_requires_reference():
    with pytest.raises(VflankError, match="ref-genome is required"):
        run_small(io.StringIO(_MAF), genome_build="hg19", ref_source="file", ref_genome=None)
