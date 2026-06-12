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

    The record is keyed on the variant identity (CHR_POS_REF_ALT), not the
    sample, so the same variant across samples deduplicates to one record:

    >{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
    {left}[REF/ALT]{right}
    >Masked__{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
    {masked_left}[REF/ALT]{masked_right}
    """
    key = f"{variant.chrom}_{variant.start}_{ref}_{alt}"
    base = "__".join(safe_header(s) for s in (variant.gene, variant.protein, variant.cdna, key))
    return [
        f">{base}\n{flanks.left}[{ref}/{alt}]{flanks.right}\n",
        f">Masked__{base}\n{flanks.masked_left}[{ref}/{alt}]{flanks.masked_right}\n",
    ]


def write_fasta(path: Path, records: list[str]) -> None:
    try:
        path.write_text("".join(records))
    except OSError as exc:
        raise OSError(f"Could not write output: {exc}") from exc
