"""Read small variants (SNP/indel) from a VCF/BCF into canonical MAF columns.

Uses ``pysam.VariantFile`` (plain / bgzip / BCF handled; no index needed for a
full sequential read). The output is a DataFrame whose columns are the *same
canonical names* as ``io/maf`` (``Chromosome``/``Start_Position``/… ), so a VCF
flows through the existing small-variant pipeline (``parse_variant_row`` →
``iter_small``) unchanged.

Conventions (see ``docs/research/sv-vcf-input.md``):

- **Sites-only.** Sample genotypes are ignored; every record gets ``SAMPLE``.
- **Multi-allelic** records expand to one row per ALT.
- **Symbolic / SV / BND** ALTs (``<DEL>``, ``A[2:123[``, ``*``, …) are skipped —
  this is the small-variant path; SV-VCF is a separate (future) feature.
- VCF's anchor-base REF/ALT encoding is normalised to MAF's 1-based, fully
  closed ``[Start, End]`` by :func:`vcf_to_maf_coords` (a pure, tested helper).
- Gene / HGVS are pulled best-effort from a VEP ``CSQ`` or SnpEff ``ANN`` INFO
  field when present, and left blank otherwise — the ``CHROM_POS_REF_ALT`` key
  in the FASTA header is always present regardless.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TypeGuard

from ..errors import VcfError
from ..logging import get_logger
from .maf import (
    MAF_ALT,
    MAF_CDNA,
    MAF_CHR,
    MAF_END,
    MAF_GENE,
    MAF_PROT,
    MAF_REF,
    MAF_SAMPLE,
    MAF_START,
)

log = get_logger()

# Extensions we treat as VCF/BCF (case-insensitive). ``.bgz`` is an alternate
# bgzip suffix some tools emit.
_VCF_SUFFIXES = (".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")

# A "simple" ALT is a plain nucleotide string; anything else (``<DEL>``, BND
# brackets, the ``*`` spanning-deletion marker) is not a small variant.
_SIMPLE_ALLELE = re.compile(r"^[ACGTNacgtn]+$")

# INFO annotation fields → the sub-keys we read gene / HGVS from. The first
# present sub-key in each tuple wins (VEP and SnpEff name things differently).
_ANNOTATORS: dict[str, dict[str, tuple[str, ...]]] = {
    "CSQ": {"gene": ("SYMBOL", "Gene"), "cdna": ("HGVSc",), "protein": ("HGVSp", "HGVSp_Short")},
    "ANN": {"gene": ("Gene_Name",), "cdna": ("HGVS.c",), "protein": ("HGVS.p",)},
}


def is_vcf_path(path: object) -> TypeGuard[str | Path]:
    """True if ``path`` looks like a VCF/BCF by file extension (case-insensitive).

    A :class:`~typing.TypeGuard` so callers narrow the input to ``str | Path`` in
    the true branch. Only filesystem paths/strings are classified; open buffers
    always read as MAF (pysam needs a real path), so they return ``False``.
    """
    if not isinstance(path, (str, Path)):
        return False
    return str(path).lower().endswith(_VCF_SUFFIXES)


def vcf_to_maf_coords(pos: int, ref: str, alt: str) -> tuple[int, int, str, str]:
    """Map a VCF ``(POS, REF, ALT)`` to MAF 1-based closed ``(start, end, ref, alt)``.

    - **SNP / MNV** (equal lengths): ``start=POS``, ``end=POS+len-1``, alleles
      unchanged.
    - **insertion** (REF is a prefix of ALT): the insertion sits *between* two
      genomic bases → ``start=anchor``, ``end=anchor+1``, ``ref='-'``,
      ``alt=<inserted bases>`` (matches the MAF insertion convention).
    - **deletion** (ALT is a prefix of REF): deleted span
      ``start=POS+len(ALT)``, ``end=POS+len(REF)-1``, ``ref=<deleted bases>``,
      ``alt='-'``.
    - **complex / non-anchored**: represented as a block substitution
      (``start=POS``, ``end=POS+len(REF)-1``, alleles unchanged).
    """
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) == len(alt):
        return pos, pos + len(ref) - 1, ref, alt
    if len(alt) > len(ref) and alt.startswith(ref):
        anchor = pos + len(ref) - 1
        return anchor, anchor + 1, "-", alt[len(ref):]
    if len(ref) > len(alt) and ref.startswith(alt):
        return pos + len(alt), pos + len(ref) - 1, ref[len(alt):], "-"
    return pos, pos + len(ref) - 1, ref, alt


def _annotation_format(header, key: str) -> list[str] | None:
    """The pipe-delimited sub-field names declared for an INFO annotation.

    VEP writes ``... Format: Allele|Consequence|…``; SnpEff writes
    ``Functional annotations: 'Allele | Annotation | …'``. Returns the field
    list, or ``None`` if the INFO key is absent / not pipe-structured.
    """
    try:
        desc = header.info[key].description
    except (KeyError, AttributeError):
        return None
    if not desc or "|" not in desc:
        return None
    _, _, spec = desc.rpartition(":")
    fields = [f.strip().strip("'\"") for f in spec.split("|")]
    return fields if len(fields) > 1 else None


def _pick(parts: list[str], idx: dict[str, int], names: tuple[str, ...]) -> str:
    """First non-empty annotation sub-field in ``names`` (by header index)."""
    for name in names:
        i = idx.get(name)
        if i is not None and i < len(parts) and parts[i]:
            return parts[i]
    return ""


def _extract_annotation(rec, alt: str, fmt_cache: dict[str, list[str] | None]) -> dict[str, str]:
    """Best-effort gene / HGVS for ``alt`` from the record's VEP/SnpEff annotation.

    Prefers the annotation entry whose leading Allele matches ``alt`` (correct
    for multi-allelic SNPs); falls back to the first entry otherwise. Empty dict
    if no annotation is present.
    """
    for key, subkeys in _ANNOTATORS.items():
        fields = fmt_cache.get(key)
        if not fields:
            continue
        raw = rec.info.get(key)
        if not raw:
            continue
        entries = list(raw) if isinstance(raw, (tuple, list)) else [raw]
        matched: list[str] | None = None
        for entry in entries:
            candidate = str(entry).split("|")
            if candidate and candidate[0] == alt:
                matched = candidate
                break
        parts: list[str] = matched if matched is not None else str(entries[0]).split("|")
        idx = {name: i for i, name in enumerate(fields)}
        return {
            "gene": _pick(parts, idx, subkeys["gene"]),
            "cdna": _pick(parts, idx, subkeys["cdna"]),
            "protein": _pick(parts, idx, subkeys["protein"]),
        }
    return {}


def load_vcf(path: str | Path, *, n_rows: int | None = None):
    """Read a VCF/BCF into a DataFrame with canonical MAF column names.

    One row per ``(record, simple ALT)``. ``n_rows`` caps the number of emitted
    rows (for previews). Raises :class:`VcfError` if the file cannot be opened
    or yields no small variants.
    """
    import pandas as pd
    import pysam

    try:
        vf = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise VcfError(f"Could not open VCF: {exc}") from exc

    fmt_cache = {k: _annotation_format(vf.header, k) for k in _ANNOTATORS}
    rows: list[dict] = []
    n_symbolic = 0
    try:
        stop = False
        for rec in vf:
            if stop:
                break
            if rec.ref is None or rec.alts is None:
                continue
            for alt in rec.alts:
                if alt is None or not _SIMPLE_ALLELE.match(alt):
                    n_symbolic += 1
                    continue
                start, end, ref_out, alt_out = vcf_to_maf_coords(rec.pos, rec.ref, alt)
                ann = _extract_annotation(rec, alt, fmt_cache)
                rows.append({
                    MAF_CHR: rec.chrom,
                    MAF_START: start,
                    MAF_END: end,
                    MAF_REF: ref_out,
                    MAF_ALT: alt_out,
                    MAF_GENE: ann.get("gene") or "UNKNOWN",
                    MAF_PROT: ann.get("protein", ""),
                    MAF_CDNA: ann.get("cdna", ""),
                    MAF_SAMPLE: "SAMPLE",
                })
                if n_rows is not None and len(rows) >= n_rows:
                    stop = True
                    break
    except Exception as exc:  # noqa: BLE001  (any read failure -> typed error)
        raise VcfError(f"Could not read VCF: {exc}") from exc
    finally:
        vf.close()

    if n_symbolic:
        log.info("Skipped %d symbolic/SV ALT allele(s) in VCF (small-variant path).", n_symbolic)
    if not rows:
        raise VcfError(
            "No small variants found in VCF (only symbolic/SV/empty records?)."
        )
    return pd.DataFrame(rows)
