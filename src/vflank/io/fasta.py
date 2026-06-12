"""FASTA record formatting and writing for the small-variant path."""

from __future__ import annotations

import re
from pathlib import Path

from ..core.flanks import FlankResult
from ..core.variant import Variant

_UNSAFE_RE = re.compile(r"[^\w.\-]")


def safe_header(s: str) -> str:
    """Replace characters unsafe in FASTA headers with underscores."""
    if not s or s in ("nan", "None", "."):
        return "."
    return _UNSAFE_RE.sub("_", s.strip())


def format_records(variant: Variant, flanks: FlankResult, ref: str, alt: str) -> list[str]:
    """Two FASTA records per variant: raw and masked.

    >{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}
    {left}[REF/ALT]{right}
    >Masked__{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}
    {masked_left}[REF/ALT]{masked_right}
    """
    base = "__".join(
        safe_header(s) for s in (variant.sample, variant.gene, variant.protein, variant.cdna)
    )
    return [
        f">{base}\n{flanks.left}[{ref}/{alt}]{flanks.right}\n",
        f">Masked__{base}\n{flanks.masked_left}[{ref}/{alt}]{flanks.masked_right}\n",
    ]


def write_fasta(path: Path, records: list[str]) -> None:
    try:
        path.write_text("".join(records))
    except OSError as exc:
        raise OSError(f"Could not write output: {exc}") from exc
