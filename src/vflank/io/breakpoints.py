"""Read SV breakpoints from the simple iCallSV / iAnnotateSV TSV.

Columns are matched **by header name, not position** (``SvColumns``), so a file
works regardless of column order or extra columns, as long as the named columns
are present. Mirrors ``io/maf.MafColumns``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..core.chrom import normalise_chrom
from ..core.fusion import Breakpoint, Fusion
from ..errors import SvError


@dataclass
class SvColumns:
    """Logical field -> header column name (all overridable)."""

    chr1: str = "chr1"
    pos1: str = "pos1"
    str1: str = "str1"
    chr2: str = "chr2"
    pos2: str = "pos2"
    str2: str = "str2"
    name: str = "name"      # optional
    sample: str = "sample"  # optional


# Logical fields that must be present (resolved through SvColumns to header names).
_REQUIRED_FIELDS = ("chr1", "pos1", "str1", "chr2", "pos2", "str2")


def read_sv_table(path: Path):
    """Read the TSV into a DataFrame (tab-separated, '#'-comment aware)."""
    import pandas as pd

    try:
        return pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    except Exception as exc:  # noqa: BLE001
        raise SvError(f"Could not read SV table: {exc}") from exc


def load_sv_table(path: Path, cols: SvColumns):
    """Read and validate that the required columns exist (by name)."""
    df = read_sv_table(path)
    if df.empty:
        raise SvError("SV table is empty.")

    missing = [getattr(cols, f) for f in _REQUIRED_FIELDS if getattr(cols, f) not in df.columns]
    if missing:
        raise SvError(
            f"Missing required column(s): {', '.join(missing)}. "
            f"Header has: {', '.join(map(str, df.columns))}. "
            "Use the column overrides to map differently-named columns."
        )
    return df


def _parse_strand(raw) -> int | None:
    try:
        value = int(float(raw))
    except (ValueError, TypeError):
        return None
    return value if value in (0, 1) else None


def _optional(row, col: str) -> str:
    import pandas as pd

    value = row.get(col)
    return str(value) if value is not None and pd.notna(value) else ""


def parse_fusion_row(row, cols: SvColumns) -> tuple[Fusion | None, str | None]:
    """Convert one row to a :class:`Fusion`, or return a skip reason."""
    c1, err1 = normalise_chrom(row[cols.chr1])
    if c1 is None:
        return None, f"breakpoint 1 — {err1}"
    c2, err2 = normalise_chrom(row[cols.chr2])
    if c2 is None:
        return None, f"breakpoint 2 — {err2}"

    try:
        p1 = int(float(row[cols.pos1]))
        p2 = int(float(row[cols.pos2]))
    except (ValueError, TypeError):
        return None, f"non-numeric position (pos1={row[cols.pos1]!r}, pos2={row[cols.pos2]!r})"
    if p1 < 1 or p2 < 1:
        return None, f"position < 1 (pos1={p1}, pos2={p2})"

    s1 = _parse_strand(row[cols.str1])
    s2 = _parse_strand(row[cols.str2])
    if s1 is None or s2 is None:
        return None, f"strand must be 0 or 1 (str1={row[cols.str1]!r}, str2={row[cols.str2]!r})"

    return (
        Fusion(
            Breakpoint(c1, p1, s1),
            Breakpoint(c2, p2, s2),
            name=_optional(row, cols.name),
            sample=_optional(row, cols.sample),
        ),
        None,
    )
