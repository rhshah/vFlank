"""Emit Primer3 Boulder-IO input from vflank's masked targets.

Primer3 reads ``TAG=value`` lines, one record per ``=``. For each variant /
junction we give it:

- ``SEQUENCE_TEMPLATE`` — the best-known sequence: the masked/consensus call at
  every position, falling back to the reference base where the call is ``N``
  (so Primer3 still has real bases for thermodynamics).
- ``SEQUENCE_TARGET`` — the variant / junction span the assay must cover.
- ``SEQUENCE_EXCLUDED_REGION`` — the masked positions (common SNPs, patient
  het/low-cov/insertion sites). A hard "do not place an oligo here" constraint,
  which is what masking means — unlike a degenerate ``N`` base, which Primer3
  may still design over.

This module is pure (plain strings in, Boulder-IO text out); the CLI collects
records during the run and calls :func:`write_primer3`. Boulder-IO has no comment
syntax, so provenance (vflank version, parameters) lives in the run report/log,
not in this file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# ddPCR-oriented starting defaults. These are *design* parameters the user will
# tune downstream; we emit a runnable minimum (short amplicons, an internal
# oligo for the probe) and nothing opinionated beyond that.
DEFAULT_GLOBAL_SETTINGS: dict[str, str] = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": "1",
    "PRIMER_PICK_INTERNAL_OLIGO": "1",  # the ddPCR probe
    "PRIMER_PICK_RIGHT_PRIMER": "1",
    "PRIMER_PRODUCT_SIZE_RANGE": "60-200",
}


@dataclass(slots=True)
class Primer3Record:
    """One Boulder-IO sequence record."""

    id: str
    template: str
    target: tuple[int, int]                       # (start, length), 0-based
    excluded: list[tuple[int, int]] = field(default_factory=list)


def resolve_template(raw: str, masked: str) -> str:
    """Best-known base per position: the masked/consensus call, or the reference
    base where the call is ``N``. ``raw`` and ``masked`` are equal length."""
    return "".join(m if m != "N" else r for r, m in zip(raw, masked, strict=False))


def n_runs(masked: str, offset: int = 0) -> list[tuple[int, int]]:
    """Coalesce runs of ``N`` in ``masked`` into ``(start + offset, length)`` pairs."""
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for i, ch in enumerate(masked):
        if ch == "N" and start is None:
            start = i
        elif ch != "N" and start is not None:
            runs.append((start + offset, i - start))
            start = None
    if start is not None:
        runs.append((start + offset, len(masked) - start))
    return runs


def small_variant_record(
    record_id: str, left: str, right: str, masked_left: str, masked_right: str, ref: str
) -> Primer3Record:
    """Build a record for a small variant. The variant (``ref`` bases) is the
    target between the two resolved flanks; masked flank positions are excluded."""
    ref_bases = ref.replace("-", "")
    tleft = resolve_template(left, masked_left)
    tright = resolve_template(right, masked_right)
    template = tleft + ref_bases + tright
    target = (len(tleft), max(1, len(ref_bases)))
    excluded = n_runs(masked_left) + n_runs(masked_right, offset=len(tleft) + len(ref_bases))
    return Primer3Record(record_id, template, target, excluded)


def junction_record(
    record_id: str, sequence: str, masked_sequence: str, junction_index: int
) -> Primer3Record:
    """Build a record for a fusion junction. The target straddles the junction so
    the amplicon spans it; masked positions are excluded."""
    template = resolve_template(sequence, masked_sequence)
    start = max(0, junction_index - 1)
    length = max(1, min(2, len(template) - start))
    return Primer3Record(record_id, template, (start, length), n_runs(masked_sequence))


def format_record(rec: Primer3Record) -> str:
    """Render one Boulder-IO record (without the trailing settings block)."""
    lines = [
        f"SEQUENCE_ID={rec.id}",
        f"SEQUENCE_TEMPLATE={rec.template}",
        f"SEQUENCE_TARGET={rec.target[0]},{rec.target[1]}",
    ]
    if rec.excluded:
        lines.append(
            "SEQUENCE_EXCLUDED_REGION=" + " ".join(f"{s},{ln}" for s, ln in rec.excluded)
        )
    return "\n".join(lines)


def write_primer3(
    path: Path, records: list[Primer3Record], settings: dict[str, str] | None = None
) -> None:
    """Write a Boulder-IO file. Global ``PRIMER_*`` settings go in the first
    record and persist for the rest (Primer3 semantics). Raises OSError on
    failure (never silent)."""
    settings = DEFAULT_GLOBAL_SETTINGS if settings is None else settings
    blocks: list[str] = []
    for i, rec in enumerate(records):
        body = format_record(rec)
        if i == 0:
            prelude = "\n".join(f"{k}={v}" for k, v in settings.items())
            body = f"{prelude}\n{body}"
        blocks.append(f"{body}\n=")
    try:
        path.write_text("\n".join(blocks) + "\n")
    except OSError as exc:
        raise OSError(f"Could not write Primer3 output: {exc}") from exc
