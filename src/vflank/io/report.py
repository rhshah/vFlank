"""Write a machine-readable TSV run report alongside the FASTA output.

Aggregate stats and the skip breakdown go in ``#``-comment header lines; the
per-variant table follows as proper TSV. Columns are taken from the row keys
(insertion order), so callers control the columns per run mode.
"""

from __future__ import annotations

from pathlib import Path


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

    columns = list(summary_rows[0].keys()) if summary_rows else []
    lines.append("\t".join(columns))
    for r in summary_rows:
        lines.append("\t".join(str(r.get(c, "")) for c in columns))

    try:
        path.write_text("\n".join(lines) + "\n")
    except OSError as exc:
        raise OSError(f"Could not write report: {exc}") from exc
