"""Write a machine-readable TSV run report alongside the FASTA output.

Aggregate stats and the skip breakdown go in ``#``-comment header lines; the
per-variant table follows as proper TSV so it loads cleanly in pandas/R.
"""

from __future__ import annotations

from pathlib import Path

_COLUMNS = ["Sample", "Gene", "Chrom", "Start", "End", "Ref", "Alt",
            "LeftLen", "RightLen", "NMasked", "Truncated"]


def write_report(
    path: Path,
    summary_rows: list[dict],
    stats: dict[str, object],
    skip_breakdown: dict[str, int],
) -> None:
    """Write the run report TSV. Raises OSError on write failure (never silent)."""
    lines: list[str] = ["# vflank small run report"]
    for key, value in stats.items():
        lines.append(f"# {key}\t{value}")
    for category, count in skip_breakdown.items():
        lines.append(f"# skip:{category}\t{count}")

    lines.append("\t".join(_COLUMNS))
    for r in summary_rows:
        lines.append("\t".join(str(r[c]) for c in _COLUMNS))

    try:
        path.write_text("\n".join(lines) + "\n")
    except OSError as exc:
        raise OSError(f"Could not write report: {exc}") from exc
