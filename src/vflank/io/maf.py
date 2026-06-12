"""MAF loading, column remapping/validation, and row -> Variant parsing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..core.chrom import normalise_chrom
from ..core.variant import Variant, validate_allele, validate_coordinates
from ..errors import MafError

# Standard TCGA/MSK MAF column names.
MAF_CHR = "Chromosome"
MAF_START = "Start_Position"
MAF_END = "End_Position"
MAF_REF = "Reference_Allele"
MAF_ALT = "Tumor_Seq_Allele2"
MAF_GENE = "Hugo_Symbol"
MAF_PROT = "HGVSp_Short"
MAF_CDNA = "HGVSc"
MAF_SAMPLE = "Tumor_Sample_Barcode"

REQUIRED_MAF_COLS = [MAF_CHR, MAF_START, MAF_END, MAF_REF, MAF_ALT]


@dataclass
class MafColumns:
    """User-overridable mapping from MAF column names to canonical names."""

    chrom: str = MAF_CHR
    start: str = MAF_START
    end: str = MAF_END
    ref: str = MAF_REF
    alt: str = MAF_ALT
    gene: str = MAF_GENE
    protein: str = MAF_PROT
    cdna: str = MAF_CDNA
    sample: str = MAF_SAMPLE


def read_maf(path: Path, *, n_rows: int | None = None):
    """Read a MAF into a DataFrame (tab-separated, '#'-comment aware)."""
    import pandas as pd

    try:
        return pd.read_csv(path, sep="\t", comment="#", nrows=n_rows, low_memory=False)
    except Exception as exc:  # noqa: BLE001
        raise MafError(f"Could not read MAF: {exc}") from exc


def load_maf(path: Path, cols: MafColumns):
    """Read a MAF, remap required columns to canonical names, and validate.

    Returns the DataFrame with canonical required-column names guaranteed
    present, and optional metadata columns filled with defaults if absent.
    """
    df = read_maf(path)
    if df.empty:
        raise MafError("MAF file is empty.")

    remap = {}
    for user_col, canonical in [
        (cols.chrom, MAF_CHR), (cols.start, MAF_START), (cols.end, MAF_END),
        (cols.ref, MAF_REF), (cols.alt, MAF_ALT),
    ]:
        if user_col != canonical and user_col in df.columns:
            remap[user_col] = canonical
    if remap:
        df = df.rename(columns=remap)

    missing = [c for c in REQUIRED_MAF_COLS if c not in df.columns]
    if missing:
        raise MafError(
            f"Missing required columns: {', '.join(missing)}. "
            "Run `vflank small inspect` to see column names, then use "
            "--chrom-col / --start-col etc. to remap."
        )

    for col, default in {
        cols.gene: "UNKNOWN", cols.protein: "", cols.cdna: "", cols.sample: "SAMPLE",
    }.items():
        if col not in df.columns:
            df[col] = default

    return df


def parse_variant_row(row, cols: MafColumns) -> tuple[Variant | None, str | None]:
    """Convert a MAF row to a :class:`Variant`, or return a skip reason.

    Returns ``(variant, None)`` on success or ``(None, reason)`` on a bad row.
    """
    import pandas as pd

    ref = str(row[MAF_REF]) if pd.notna(row[MAF_REF]) else "-"
    alt = str(row[MAF_ALT]) if pd.notna(row[MAF_ALT]) else "-"
    gene = str(row.get(cols.gene, "UNKNOWN"))

    bare, chrom_err = normalise_chrom(row[MAF_CHR])
    if bare is None:
        # By contract, normalise_chrom sets chrom_err whenever bare is None.
        return None, f"{gene} — {chrom_err}"

    raw = str(row[MAF_CHR]).strip()
    raw_display = raw if raw.upper().startswith("CHR") else bare

    try:
        start = int(float(row[MAF_START]))
        end = int(float(row[MAF_END]))
    except (ValueError, TypeError):
        return None, (
            f"{gene}: non-numeric position "
            f"(start={row[MAF_START]!r}, end={row[MAF_END]!r})"
        )

    coord_err = validate_coordinates(start, end)
    if coord_err:
        return None, f"{gene} {raw_display}:{start} — {coord_err}"

    if not validate_allele(ref) or not validate_allele(alt):
        return None, f"{gene} {raw_display}:{start} — invalid allele (ref={ref!r}, alt={alt!r})"

    return (
        Variant(
            chrom=bare, start=start, end=end, ref=ref, alt=alt,
            gene=gene,
            protein=str(row.get(cols.protein, "")),
            cdna=str(row.get(cols.cdna, "")),
            sample=str(row.get(cols.sample, "SAMPLE")),
            raw_chrom=raw_display,
        ),
        None,
    )
