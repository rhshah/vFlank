import pytest

from vflank.errors import SvError
from vflank.io.breakpoints import (
    SvColumns,
    load_sv_table,
    parse_fusion_row,
)


def _write(tmp_path, header, rows):
    p = tmp_path / "sv.txt"
    lines = ["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]
    p.write_text("\n".join(lines) + "\n")
    return p


def test_columns_matched_by_name_any_order(tmp_path):
    # Columns in a scrambled order, with an extra unrelated column.
    header = ["sample", "pos2", "str2", "chr1", "extra", "pos1", "str1", "chr2", "name"]
    rows = [["S1", 32416521, 1, 22, "junk", 29684066, 0, 11, "EWSR1-WT1"]]
    df = load_sv_table(_write(tmp_path, header, rows), SvColumns())
    fus, reason = parse_fusion_row(df.iloc[0], SvColumns())
    assert reason is None
    assert fus.bp1.chrom == "22" and fus.bp1.pos == 29684066 and fus.bp1.strand == 0
    assert fus.bp2.chrom == "11" and fus.bp2.pos == 32416521 and fus.bp2.strand == 1
    assert fus.name == "EWSR1-WT1" and fus.sample == "S1"


def test_renamed_columns_via_override(tmp_path):
    header = ["chromosome1", "position1", "strand1", "chromosome2", "position2", "strand2"]
    rows = [[7, 100, 0, 7, 200, 0]]
    cols = SvColumns(
        chr1="chromosome1", pos1="position1", str1="strand1",
        chr2="chromosome2", pos2="position2", str2="strand2",
    )
    df = load_sv_table(_write(tmp_path, header, rows), cols)
    fus, reason = parse_fusion_row(df.iloc[0], cols)
    assert reason is None and fus.bp1.pos == 100 and fus.bp2.pos == 200


def test_missing_required_column_errors(tmp_path):
    header = ["chr1", "pos1", "str1", "chr2", "pos2"]  # no str2
    with pytest.raises(SvError, match="str2"):
        load_sv_table(_write(tmp_path, header, [[1, 10, 0, 2, 20]]), SvColumns())


def test_bad_strand_and_position_are_skipped(tmp_path):
    cols = SvColumns()
    header = ["chr1", "pos1", "str1", "chr2", "pos2", "str2"]
    df = load_sv_table(
        _write(tmp_path, header, [
            [1, 10, 2, 2, 20, 0],     # str1=2 invalid
            [1, -5, 0, 2, 20, 1],     # pos1 < 1
        ]),
        cols,
    )
    _, r0 = parse_fusion_row(df.iloc[0], cols)
    _, r1 = parse_fusion_row(df.iloc[1], cols)
    assert r0 and "strand" in r0
    assert r1 and "position" in r1
