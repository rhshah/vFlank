"""Unit tests for the presentation-free orchestration (vflank.pipeline).

These exercise iter_small / collect directly with a tiny fake reference — no
pysam, no CliRunner — which is the whole point of the extraction: the
orchestration is testable without the CLI.
"""

from __future__ import annotations

import pandas as pd

from vflank.io.maf import MafColumns
from vflank.pipeline import Duplicate, Processed, Skipped, collect, iter_small

SEQ = "".join("ACGT"[i % 4] for i in range(60))


class FakeReference:
    """Minimal reference: fetch(chrom, start0, end0) -> slice (0-based half-open,
    pysam-style; a window past the end yields a short string, like a contig edge)."""

    has_chr = False

    def __init__(self, seq: str) -> None:
        self.seq = seq

    def fetch(self, chrom: str, s0: int, e0: int) -> str:
        return self.seq[s0:e0]

    def close(self) -> None:
        pass


_MAF_COLS = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
             "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]


def _df(rows: list[list]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=_MAF_COLS)


def _row(gene="TP53", chrom="1", start=30, end=30, ref="A", alt="T", sample="S1"):
    return [gene, chrom, start, end, ref, alt, sample]


def _iter(df, **kw):
    opts = dict(cols=MafColumns(), reference=FakeReference(SEQ), gnomad=None,
                flank=5, af_threshold=0.001, dedup=True, uppercase=True)
    opts.update(kw)
    return iter_small(df, **opts)


def test_processed_emits_records_and_detail():
    outcomes = list(_iter(_df([_row()])))
    assert len(outcomes) == 1
    o = outcomes[0]
    assert isinstance(o, Processed)
    assert o.records[0].startswith(">TP53__")
    assert o.detail["LeftLen"] == 5 and o.detail["RightLen"] == 5
    assert o.truncation is None
    assert o.primer3 is None  # emit_primer3 defaults off


def test_bad_row_is_skipped_not_dropped():
    outcomes = list(_iter(_df([_row(start="bad")])))  # non-numeric position
    assert len(outcomes) == 1 and isinstance(outcomes[0], Skipped)
    assert outcomes[0].message.startswith("row 0")


def test_dedup_yields_duplicate():
    outcomes = list(_iter(_df([_row(), _row()])))  # same variant twice
    assert [type(o).__name__ for o in outcomes] == ["Processed", "Duplicate"]
    assert isinstance(outcomes[1], Duplicate)


def test_truncation_flagged_but_emitted():
    # Variant at 1-based pos 58: right flank (58..63) runs off the 60 bp contig.
    o = next(iter(_iter(_df([_row(start=58, end=58)]))))
    assert isinstance(o, Processed)                 # still emitted
    assert o.truncation is not None and "truncated" in o.truncation
    assert o.detail["Truncated"] is True


def test_emit_primer3_builds_record():
    o = next(iter(_iter(_df([_row()]), emit_primer3=True)))
    assert isinstance(o, Processed) and o.primer3 is not None
    assert o.primer3.id.startswith("TP53__")


def test_collect_accumulates_counts_and_breakdown():
    df = _df([_row(), _row(), _row(start="bad")])
    result = collect(_iter(df))
    assert result.n_processed == 1
    assert result.n_duplicate == 1
    assert result.n_skipped == 1
    assert sum(result.skip_breakdown.values()) == 1
    assert result.skip_breakdown.most_common(1)  # it's a Counter (CLI summary uses this)
    assert len(result.records) == 2  # raw + masked for the one processed variant
