from vflank.io.report import write_report


def test_write_report(tmp_path):
    rows = [{
        "Sample": "S1", "Gene": "TP53", "Chrom": "17", "Start": 100, "End": 100,
        "Ref": "A", "Alt": "T", "LeftLen": 200, "RightLen": 200,
        "NMasked": 3, "Truncated": False,
    }]
    stats = {"total_in_maf": 2, "processed": 1, "skipped": 1}
    out = tmp_path / "report.tsv"
    write_report(out, rows, stats, {"chromosome missing/invalid": 1})

    text = out.read_text()
    lines = text.splitlines()
    assert "# total_in_maf\t2" in lines
    assert "# skip:chromosome missing/invalid\t1" in lines
    # Header line then one data row, both proper TSV.
    header = next(line for line in lines if line.startswith("Sample\t"))
    assert header.split("\t") == [
        "Sample", "Gene", "Chrom", "Start", "End", "Ref", "Alt",
        "LeftLen", "RightLen", "NMasked", "Truncated",
    ]
    data = lines[-1].split("\t")
    assert data[0] == "S1" and data[3] == "100" and data[9] == "3"
