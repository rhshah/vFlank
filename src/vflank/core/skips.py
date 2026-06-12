"""Categorise per-variant skip reasons for a compact run summary.

The reasons are free-text (they come from several validation points), so this
buckets them by keyword into a small, stable set of categories. Pure and
testable.
"""

from __future__ import annotations

# Order matters: the first matching keyword wins.
_CATEGORIES: list[tuple[str, str]] = [
    ("chromosome", "chromosome missing/invalid"),
    ("non-numeric position", "non-numeric position"),
    ("invalid allele", "invalid allele"),
    ("fetch error", "fetch error"),
    ("< 1", "invalid coordinates"),
    ("< start", "invalid coordinates"),
]


def categorize_skip(reason: str) -> str:
    """Map a free-text skip reason to a stable category label."""
    r = reason.lower()
    for keyword, label in _CATEGORIES:
        if keyword in r:
            return label
    return "other"
