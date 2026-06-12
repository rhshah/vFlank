"""
get_flanking_sequence.py
~~~~~~~~~~~~~~~~~~~~~~~~

Given a MAF file, retrieves N bp of flanking sequence around each variant from
a reference genome FASTA, masks positions that are common germline SNPs in a
population VCF (gnomAD), and writes a FASTA file with both raw and masked
records per variant — suitable for ddPCR assay design.

Supports hg19 (GRCh37) and hg38 (GRCh38).

Chromosome notation is auto-detected from the FASTA and tabix indices at
startup, so the tool works regardless of whether the files use 'chr1' or '1'
notation, and regardless of what the MAF contains.

Population VCF directory layout
---------------------------------
gnomAD v4.1 (GRCh38, per-chromosome — recommended):
  Directory should contain files named:
    gnomad.genomes.v4.1.sites.chr1.vcf.bgz  (+ .tbi)
    ...
    gnomad.genomes.v4.1.sites.chrX.vcf.bgz
  Download:
    for CHR in {1..22} X Y; do
      wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/
           gnomad.genomes.v4.1.sites.chr${CHR}.vcf.bgz{,.tbi}
    done

gnomAD v2.1.1 (GRCh37/hg19, per-chromosome):
  Directory should contain files named:
    gnomad.genomes.r2.1.1.sites.1.vcf.bgz  (+ .tbi)
    ...
  Download:
    for CHR in {1..22} X Y; do
      wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/
           gnomad.genomes.r2.1.1.sites.${CHR}.vcf.bgz{,.tbi}
    done

Subcommands
-----------
  run      – main analysis
  inspect  – preview MAF columns / first rows
  list-vcf – verify per-chromosome VCF coverage in a directory

Example
-------
  # All samples
  python get_flanking_sequence.py run variants.maf \\
      --ref-genome /path/to/GRCh38.fasta \\
      --pop-vcf-dir /path/to/gnomad_v4/ \\
      --genome-build hg38 \\
      --output flanking_sequences.fasta

  # One sample by ID
  python get_flanking_sequence.py run variants.maf \\
      --ref-genome /path/to/GRCh38.fasta \\
      --samples "SAMPLE_001" ...

  # Multiple samples, comma-separated
  python get_flanking_sequence.py run variants.maf \\
      --samples "SAMPLE_001,SAMPLE_002,SAMPLE_003" ...

  # Sample list from file (one ID per line, # comments allowed)
  python get_flanking_sequence.py run variants.maf \\
      --samples-file /path/to/sample_list.txt ...
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Set

import typer
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich import print as rprint

# ---------------------------------------------------------------------------
# Lazy imports
# ---------------------------------------------------------------------------

def _require_pysam():
    try:
        import pysam
        return pysam
    except ImportError:
        rprint("[bold red]ERROR:[/bold red] pysam is not installed.\n"
               "  Install with: [cyan]pip install pysam[/cyan]")
        raise typer.Exit(1)


def _require_pandas():
    try:
        import pandas as pd
        return pd
    except ImportError:
        rprint("[bold red]ERROR:[/bold red] pandas is not installed.\n"
               "  Install with: [cyan]pip install pandas[/cyan]")
        raise typer.Exit(1)


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = typer.Typer(
    name="get_flanking_sequence",
    help="Extract and mask flanking sequences from a MAF file for ddPCR assay design.",
    add_completion=False,
    no_args_is_help=True,
)
console = Console(stderr=True)

# ---------------------------------------------------------------------------
# MAF column constants (standard TCGA/MSK MAF)
# ---------------------------------------------------------------------------

MAF_CHR    = "Chromosome"
MAF_START  = "Start_Position"
MAF_END    = "End_Position"
MAF_REF    = "Reference_Allele"
MAF_ALT    = "Tumor_Seq_Allele2"
MAF_GENE   = "Hugo_Symbol"
MAF_PROT   = "HGVSp_Short"
MAF_CDNA   = "HGVSc"
MAF_SAMPLE = "Tumor_Sample_Barcode"

REQUIRED_MAF_COLS = [MAF_CHR, MAF_START, MAF_END, MAF_REF, MAF_ALT]

VALID_BASES  = set("ATGCatgc-N")
VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

# gnomAD filename patterns, keyed by genome build.
# {chrom} is substituted with the bare chromosome (e.g. "1", "X") at runtime.
GNOMAD_PATTERNS: Dict[str, List[str]] = {
    "hg38": [
        "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz",
        "gnomad.genomes.v4.0.sites.chr{chrom}.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr{chrom}.vcf.bgz",
        "gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz",
    ],
    "hg19": [
        "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz",
        "gnomad.genomes.r2.1.2.sites.{chrom}.vcf.bgz",
        "ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz",   # single-file legacy fallback
    ],
}

# ---------------------------------------------------------------------------
# Chromosome notation detection and normalisation
# ---------------------------------------------------------------------------

# Numeric encodings some tools use for sex / mito chromosomes
_NUMERIC_CHROM_MAP: Dict[str, str] = {
    "23": "X",
    "24": "Y",
    "25": "MT",
    "26": "MT",
}

# Common mitochondrial aliases → canonical "MT"
_MITO_ALIASES: Dict[str, str] = {
    "M":   "MT",
    "MT":  "MT",
    "chrM": "MT",
    "chrMT": "MT",
}


def _detect_fasta_chr_style(fa) -> bool:
    """
    Return True if the FASTA uses 'chr'-prefixed contig names (e.g. 'chr1'),
    False if it uses bare names (e.g. '1').
    Probes chr1/1 first, then falls back through chr2-chr5/2-5.
    """
    refs = set(fa.references)
    if "chr1" in refs:
        return True
    if "1" in refs:
        return False
    for i in range(2, 6):
        if f"chr{i}" in refs:
            return True
        if str(i) in refs:
            return False
    return True   # can't determine — assume prefixed


def _detect_tabix_chr_style(tbx) -> bool:
    """Return True if the tabix VCF uses 'chr'-prefixed contig names."""
    return any(c.startswith("chr") for c in tbx.contigs)


def _detect_maf_chr_style(df, chrom_col: str) -> Optional[bool]:
    """
    Inspect the MAF Chromosome column and return:
      True   – MAF uses chr-prefixed notation (e.g. 'chr7')
      False  – MAF uses bare notation (e.g. '7')
      None   – column absent or all values are NaN / unrecognisable

    This is used only for reporting; the actual per-row normalisation is
    handled by _normalise_maf_chrom() which is unambiguous.
    """
    if chrom_col not in df.columns:
        return None
    sample = df[chrom_col].dropna().astype(str).head(20)
    for val in sample:
        v = val.strip()
        if v.lower().startswith("chr"):
            return True
        if v in VALID_CHROMS or v in _NUMERIC_CHROM_MAP:
            return False
    return None


def _normalise_maf_chrom(raw_chrom) -> tuple:
    """
    Convert a raw MAF Chromosome value to a canonical bare chrom string.

    Handles:
      - NaN / None / empty → error
      - Any case of 'chr' prefix: 'Chr7', 'CHR7', 'chr7' → '7'
      - Bare integers: '7' → '7'
      - Numeric sex/mito encodings: '23'→'X', '24'→'Y', '25'/'26'→'MT'
      - Mito aliases: 'M', 'chrM', 'chrMT' → 'MT'
      - Validation against VALID_CHROMS

    Returns: (bare_chrom: str, error_msg: str | None)
      If error_msg is not None, bare_chrom is None and the variant should be skipped.
    """
    # Guard against NaN cells (pandas reads missing values as float NaN)
    if raw_chrom is None:
        return None, "chromosome is None"
    try:
        import math
        if isinstance(raw_chrom, float) and math.isnan(raw_chrom):
            return None, "chromosome is NaN (missing value in MAF)"
    except (TypeError, ValueError):
        pass

    chrom = str(raw_chrom).strip()

    if not chrom or chrom.lower() in ("nan", "none", ".", ""):
        return None, f"chromosome is empty or missing (got {raw_chrom!r})"

    # Mitochondrial aliases first (before chr-stripping, to catch 'chrM')
    upper = chrom.upper()
    for alias, canonical in _MITO_ALIASES.items():
        if upper == alias.upper():
            return canonical, None

    # Case-insensitive chr-prefix stripping
    if upper.startswith("CHR"):
        bare = chrom[3:]   # preserve original case of the rest (e.g. 'X', 'Y')
    else:
        bare = chrom

    # Numeric alternative encodings (23→X, 24→Y, 25/26→MT)
    if bare in _NUMERIC_CHROM_MAP:
        bare = _NUMERIC_CHROM_MAP[bare]

    # Uppercase for X/Y consistency
    bare = bare.upper() if bare.upper() in ("X", "Y", "MT") else bare

    if bare not in VALID_CHROMS:
        return None, (
            f"unrecognised chromosome value {raw_chrom!r} "
            f"(normalised to {bare!r}, not in {sorted(VALID_CHROMS)})"
        )

    return bare, None


def _chrom_for_fasta(bare: str, fasta_has_chr: bool) -> str:
    """Convert bare chromosome to the notation used by the FASTA."""
    return f"chr{bare}" if fasta_has_chr else bare


def _chrom_for_tabix(bare: str, tabix_has_chr: bool) -> str:
    """Convert bare chromosome to the notation used by the tabix VCF."""
    return f"chr{bare}" if tabix_has_chr else bare


def _bare_chrom(chrom: str) -> str:
    """
    Strip leading 'chr' (case-insensitive) to get a bare chromosome string.
    Use _normalise_maf_chrom() for MAF values; this is for already-trusted strings.
    """
    if chrom.upper().startswith("CHR"):
        return chrom[3:]
    return chrom


# ---------------------------------------------------------------------------
# Per-chromosome VCF resolver
# ---------------------------------------------------------------------------

def _resolve_vcf_for_chrom(
    pop_vcf_dir: Path,
    bare_chrom: str,
    genome_build: str,
) -> Optional[Path]:
    """
    Scan GNOMAD_PATTERNS for the given build and return the first matching
    file path that exists in pop_vcf_dir.  Returns None if nothing matches.
    """
    for pattern in GNOMAD_PATTERNS.get(genome_build, []):
        if "{chrom}" not in pattern:
            # Single-file fallback (e.g. ExAC) — covers all chromosomes
            candidate = pop_vcf_dir / pattern
            if candidate.exists():
                return candidate
        else:
            candidate = pop_vcf_dir / pattern.format(chrom=bare_chrom)
            if candidate.exists():
                return candidate
    return None


def _build_chrom_vcf_map(
    pop_vcf_dir: Path,
    genome_build: str,
    chroms: List[str],
) -> Dict[str, Optional[Path]]:
    return {c: _resolve_vcf_for_chrom(pop_vcf_dir, c, genome_build) for c in chroms}


# ---------------------------------------------------------------------------
# Allele / coordinate validation
# ---------------------------------------------------------------------------

def _validate_allele(allele: str) -> bool:
    if not allele or allele in (".", "-"):
        return True
    return all(b in VALID_BASES for b in allele)


def _validate_coordinates(start: int, end: int) -> Optional[str]:
    if start < 1:
        return f"start={start} is < 1"
    if end < start:
        return f"end={end} < start={start}"
    return None


# ---------------------------------------------------------------------------
# Sample filter loader
# ---------------------------------------------------------------------------

def _load_sample_filter(
    samples_str: Optional[str],
    samples_file: Optional[Path],
) -> Optional[Set[str]]:
    """
    Parse sample IDs from --samples (comma-separated) and/or --samples-file
    (one ID per line, # comments allowed).

    Returns a set of IDs to keep, or None if no filtering is requested.
    Raises typer.Exit(1) on any file-read error.
    """
    ids: Set[str] = set()

    if samples_str:
        for token in samples_str.split(","):
            s = token.strip()
            if s:
                ids.add(s)

    if samples_file is not None:
        if not samples_file.exists():
            rprint(f"[bold red]ERROR:[/bold red] --samples-file not found: {samples_file}")
            raise typer.Exit(1)
        try:
            for line in samples_file.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    ids.add(line)
        except OSError as exc:
            rprint(f"[bold red]ERROR:[/bold red] Could not read --samples-file: {exc}")
            raise typer.Exit(1)

    if not ids and (samples_str is not None or samples_file is not None):
        rprint("[bold red]ERROR:[/bold red] --samples / --samples-file produced an empty ID set.")
        raise typer.Exit(1)

    return ids if ids else None


# ---------------------------------------------------------------------------
# Population VCF querying
# ---------------------------------------------------------------------------

def _get_pop_snp_positions(
    tabix,
    chrom_tabix: str,
    region_start_0based: int,
    region_end_0based: int,
    af_threshold: float,
) -> List[int]:
    """
    Query a tabix-indexed gnomAD VCF and return 1-based positions of SNPs
    whose maximum population AF (AF or AF_grpmax) >= af_threshold.

    Only single-base substitutions are considered for masking.
    """
    positions: List[int] = []

    try:
        rows = list(tabix.fetch(chrom_tabix, region_start_0based, region_end_0based))
    except (ValueError, KeyError):
        # Contig absent in this VCF (e.g. MT in some builds) — skip silently
        return positions

    for row in rows:
        fields = row.split("\t")
        if len(fields) < 8:
            continue

        pos       = int(fields[1])   # 1-based VCF position
        ref       = fields[3]
        alt_field = fields[4]
        info      = fields[7]

        # SNPs only: single-base REF, all ALT alleles single-base
        if len(ref) != 1:
            continue
        alts = alt_field.split(",")
        if not all(len(a) == 1 and a not in (".", "*") for a in alts):
            continue

        # Parse AF and AF_grpmax, guarding against "." and "NaN"
        af_values: List[float] = []
        for token in info.split(";"):
            if token.startswith("AF=") or token.startswith("AF_grpmax="):
                raw = token.split("=", 1)[1]
                for v in raw.split(","):
                    try:
                        f = float(v)
                        if f == f:   # reject NaN (NaN != NaN)
                            af_values.append(f)
                    except ValueError:
                        pass

        if af_values and max(af_values) >= af_threshold:
            positions.append(pos)

    return positions


# ---------------------------------------------------------------------------
# Sequence masking
# ---------------------------------------------------------------------------

def _mask_sequence(
    seq: str,
    region_start_0based: int,
    positions_1based: List[int],
) -> str:
    """
    Replace bases in seq with 'N' at the given 1-based genomic positions.
    region_start_0based is the 0-based start passed to pysam.fetch(),
    so seq[0] corresponds to genomic position (region_start_0based + 1).
    """
    if not positions_1based:
        return seq
    seq_list = list(seq)
    for pos in positions_1based:
        idx = pos - region_start_0based - 1   # 1-based → 0-based index in seq
        if 0 <= idx < len(seq_list):
            seq_list[idx] = "N"
    return "".join(seq_list)


# ---------------------------------------------------------------------------
# Core flank fetcher
# ---------------------------------------------------------------------------

def _fetch_flanks(
    fa,
    tabix_cache: Dict[str, Optional[object]],
    fasta_has_chr: bool,
    tabix_has_chr_cache: Dict[str, bool],
    bare_chrom: str,
    start_1based: int,
    end_1based: int,
    flank: int,
    af_threshold: float,
) -> tuple:
    """
    Fetch left and right flanks from the reference FASTA, then mask common
    SNPs from the per-chromosome tabix VCF.

    MAF coordinates are 1-based fully-closed [start, end].
    pysam uses 0-based half-open [start, end).

    Returns: (flank_left, flank_right, masked_left, masked_right)
    """
    chrom_fa = _chrom_for_fasta(bare_chrom, fasta_has_chr)

    # 0-based half-open intervals
    left_start_0  = max(0, start_1based - flank - 1)
    left_end_0    = start_1based - 1      # bases up to (not including) variant
    right_start_0 = end_1based            # bases starting immediately after variant
    right_end_0   = end_1based + flank

    flank_left  = fa.fetch(chrom_fa, left_start_0, left_end_0)
    flank_right = fa.fetch(chrom_fa, right_start_0, right_end_0)

    tabix = tabix_cache.get(bare_chrom)
    if tabix is None:
        return flank_left, flank_right, flank_left, flank_right

    tabix_has_chr = tabix_has_chr_cache.get(bare_chrom, False)
    chrom_tbx = _chrom_for_tabix(bare_chrom, tabix_has_chr)

    left_snps  = _get_pop_snp_positions(tabix, chrom_tbx, left_start_0,  left_end_0,  af_threshold)
    right_snps = _get_pop_snp_positions(tabix, chrom_tbx, right_start_0, right_end_0, af_threshold)

    masked_left  = _mask_sequence(flank_left,  left_start_0,  left_snps)
    masked_right = _mask_sequence(flank_right, right_start_0, right_snps)

    return flank_left, flank_right, masked_left, masked_right


# ---------------------------------------------------------------------------
# FASTA header sanitiser
# ---------------------------------------------------------------------------

_UNSAFE_RE = re.compile(r"[^\w.\-]")

def _safe(s: str) -> str:
    """Replace characters that are unsafe in FASTA headers with underscores."""
    if not s or s in ("nan", "None", "."):
        return "."
    return _UNSAFE_RE.sub("_", s.strip())


# ---------------------------------------------------------------------------
# Subcommand: inspect
# ---------------------------------------------------------------------------

@app.command("inspect")
def inspect_maf(
    maf_file: Path = typer.Argument(..., help="MAF file to preview", exists=True),
    n_rows: int = typer.Option(3, "--rows", "-n", help="Number of data rows to show", min=1, max=20),
):
    """
    Print column names and a data preview without running any analysis.
    Useful to verify column naming before a run.
    """
    pd = _require_pandas()
    df = pd.read_csv(maf_file, sep="\t", comment="#", nrows=n_rows, low_memory=False)

    console.print(Panel(
        f"[bold cyan]{maf_file.name}[/bold cyan]  "
        f"[dim]{len(df.columns)} columns · preview of {len(df)} rows[/dim]",
        title="MAF Inspect",
    ))

    table = Table(show_header=True, header_style="bold magenta", show_lines=True)
    table.add_column("#", style="dim", width=4)
    table.add_column("Column")
    table.add_column("Required?", justify="center", width=10)
    for i in range(min(n_rows, len(df))):
        table.add_column(f"Row {i+1}", overflow="fold")

    for i, col in enumerate(df.columns, 1):
        req = "[green]✓[/green]" if col in REQUIRED_MAF_COLS else ""
        row_vals = [str(df[col].iloc[j]) if j < len(df) else "" for j in range(min(n_rows, len(df)))]
        table.add_row(str(i), col, req, *row_vals)

    console.print(table)

    missing = [c for c in REQUIRED_MAF_COLS if c not in df.columns]
    if missing:
        console.print(f"\n[bold red]Missing required columns:[/bold red] {', '.join(missing)}")
    else:
        console.print("\n[bold green]✓ All required columns present.[/bold green]")

    # Chromosome notation report
    maf_chr_style = _detect_maf_chr_style(df, MAF_CHR)
    if maf_chr_style is True:
        console.print("\n[bold]Chromosome notation:[/bold] chr-prefixed (e.g. chr7)")
    elif maf_chr_style is False:
        console.print("\n[bold]Chromosome notation:[/bold] bare (e.g. 7)")
    else:
        console.print("\n[bold]Chromosome notation:[/bold] [yellow]unknown / mixed — check column[/yellow]")

    if MAF_SAMPLE in df.columns:
        samples = df[MAF_SAMPLE].dropna().unique().tolist()
        console.print(
            f"\n[bold]Samples found ({len(samples)}):[/bold] "
            + ", ".join(str(s) for s in samples[:10])
            + (" …" if len(samples) > 10 else "")
        )


# ---------------------------------------------------------------------------
# Subcommand: list-vcf
# ---------------------------------------------------------------------------

@app.command("list-vcf")
def list_vcf(
    pop_vcf_dir: Path = typer.Argument(
        ..., help="Directory containing gnomAD per-chromosome VCF bgz files", exists=True,
    ),
    genome_build: str = typer.Option("hg38", "--genome-build", "-g", help="hg38 or hg19"),
):
    """
    Show which per-chromosome VCFs were found (and which are missing).
    Run this before a full analysis to verify your gnomAD directory.
    """
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    vcf_map = _build_chrom_vcf_map(pop_vcf_dir, genome_build, chroms)

    table = Table(
        show_header=True, header_style="bold cyan",
        title=f"VCF files in {pop_vcf_dir}  [dim]({genome_build})[/dim]",
    )
    table.add_column("Chrom", width=8)
    table.add_column("File found", overflow="fold")
    table.add_column("TBI index", justify="center", width=12)

    found = 0
    for chrom in chroms:
        vcf_path = vcf_map.get(chrom)
        if vcf_path:
            tbi = Path(str(vcf_path) + ".tbi")
            tbi_str = "[green]✓[/green]" if tbi.exists() else "[red]✗ missing[/red]"
            table.add_row(chrom, vcf_path.name, tbi_str)
            found += 1
        else:
            table.add_row(chrom, "[dim]not found[/dim]", "[dim]—[/dim]")

    console.print(table)
    console.print(f"\n{found}/{len(chroms)} chromosomes have a matching VCF.")
    if found == 0:
        console.print(
            "[yellow]⚠  No files matched known gnomAD filename patterns.\n"
            "   Expected patterns for hg38 (example):\n"
            "     gnomad.genomes.v4.1.sites.chr1.vcf.bgz\n"
            "   Expected patterns for hg19 (example):\n"
            "     gnomad.genomes.r2.1.1.sites.1.vcf.bgz[/yellow]"
        )


# ---------------------------------------------------------------------------
# Subcommand: run
# ---------------------------------------------------------------------------

@app.command("run")
def run(
    maf_file: Path = typer.Argument(
        ..., help="Input MAF file (tab-separated, TCGA/MSK format)", exists=True,
    ),
    ref_genome: Path = typer.Option(
        ..., "--ref-genome", "-r",
        help="Path to indexed reference FASTA (.fai must exist alongside it)",
    ),
    pop_vcf_dir: Optional[Path] = typer.Option(
        None, "--pop-vcf-dir", "-d",
        help=(
            "Directory of per-chromosome gnomAD VCF bgz files. "
            "The correct file for each chromosome is resolved automatically "
            "from known gnomAD filename patterns. "
            "Run [list-vcf] to verify coverage. Omit to skip masking."
        ),
    ),
    genome_build: str = typer.Option(
        "hg38", "--genome-build", "-g",
        help="Reference genome build: hg38 (GRCh38) or hg19 (GRCh37)",
        show_default=True,
    ),
    flank: int = typer.Option(
        200, "--flank", "-f",
        help="Bases to include on each side of the variant",
        show_default=True, min=1, max=10_000,
    ),
    af_threshold: float = typer.Option(
        0.001, "--af-threshold",
        help="Min population AF (or AF_grpmax) to mask a SNP [default: 0.1%%]",
        show_default=True, min=0.0, max=1.0,
    ),
    output: Path = typer.Option(
        Path("flanking_sequences.fasta"), "--output", "-o",
        help="Output FASTA file",
        show_default=True,
    ),
    # --- Sample filtering ---
    samples: Optional[str] = typer.Option(
        None, "--samples", "-s",
        help=(
            "Comma-separated list of Tumor_Sample_Barcode IDs to include. "
            "Example: --samples 'P-001,P-002,P-003'"
        ),
    ),
    samples_file: Optional[Path] = typer.Option(
        None, "--samples-file",
        help=(
            "Path to a plain-text file with one sample ID per line. "
            "Lines starting with '#' are treated as comments and ignored. "
            "Can be combined with --samples; the union is used."
        ),
    ),
    # --- Column overrides ---
    chrom_col:  str = typer.Option(MAF_CHR,    "--chrom-col"),
    start_col:  str = typer.Option(MAF_START,  "--start-col"),
    end_col:    str = typer.Option(MAF_END,    "--end-col"),
    ref_col:    str = typer.Option(MAF_REF,    "--ref-col"),
    alt_col:    str = typer.Option(MAF_ALT,    "--alt-col"),
    gene_col:   str = typer.Option(MAF_GENE,   "--gene-col"),
    prot_col:   str = typer.Option(MAF_PROT,   "--prot-col"),
    cdna_col:   str = typer.Option(MAF_CDNA,   "--cdna-col"),
    sample_col: str = typer.Option(MAF_SAMPLE, "--sample-col"),
    # --- Behaviour ---
    skip_invalid: bool = typer.Option(
        True, "--skip-invalid/--fail-invalid",
        help="Skip rows with bad alleles/coords rather than aborting",
    ),
    uppercase: bool = typer.Option(
        True, "--uppercase/--no-uppercase",
        help="Convert flanking sequences to uppercase",
        show_default=True,
    ),
):
    """
    Extract flanking sequences for every variant in a MAF file and write a FASTA.

    Each variant produces two FASTA records:

    \b
      >{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}
      {left_flank}[REF/ALT]{right_flank}

      >Masked__{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}
      {left_flank_masked}[REF/ALT]{right_flank_masked}

    Common population SNPs (AF >= --af-threshold) are replaced with 'N' in
    the Masked record. Chromosome notation ('chr1' vs '1') is auto-detected
    from the FASTA and tabix indices — no manual adjustment needed.
    """
    pysam = _require_pysam()
    pd    = _require_pandas()

    t0 = time.time()
    console.rule("[bold blue]get_flanking_sequence[/bold blue]")

    # ------------------------------------------------------------------
    # Validate static inputs
    # ------------------------------------------------------------------
    if genome_build not in ("hg19", "hg38"):
        rprint(f"[bold red]ERROR:[/bold red] --genome-build must be 'hg19' or 'hg38', got '{genome_build}'")
        raise typer.Exit(1)

    if not ref_genome.exists():
        rprint(f"[bold red]ERROR:[/bold red] Reference FASTA not found: {ref_genome}")
        raise typer.Exit(1)
    fai = Path(str(ref_genome) + ".fai")
    if not fai.exists():
        rprint(
            f"[bold red]ERROR:[/bold red] FASTA index not found: {fai}\n"
            f"  Fix: [cyan]samtools faidx {ref_genome}[/cyan]"
        )
        raise typer.Exit(1)

    if pop_vcf_dir is not None and not pop_vcf_dir.is_dir():
        rprint(f"[bold red]ERROR:[/bold red] --pop-vcf-dir is not a directory: {pop_vcf_dir}")
        raise typer.Exit(1)

    # ------------------------------------------------------------------
    # Sample filter
    # ------------------------------------------------------------------
    sample_filter: Optional[Set[str]] = _load_sample_filter(samples, samples_file)
    if sample_filter is not None:
        console.print(
            f"[bold]Sample filter:[/bold] {len(sample_filter)} ID(s) — "
            + ", ".join(sorted(sample_filter)[:6])
            + (" …" if len(sample_filter) > 6 else "")
        )

    # ------------------------------------------------------------------
    # Load and filter MAF
    # ------------------------------------------------------------------
    console.print(f"[bold]Loading MAF:[/bold] {maf_file}")
    try:
        df = pd.read_csv(maf_file, sep="\t", comment="#", low_memory=False)
    except Exception as exc:
        rprint(f"[bold red]ERROR:[/bold red] Could not read MAF: {exc}")
        raise typer.Exit(1)

    if df.empty:
        rprint("[bold red]ERROR:[/bold red] MAF file is empty.")
        raise typer.Exit(1)

    # Remap user-specified column names to internal canonical names
    col_remap = {}
    for user_col, canonical in [
        (chrom_col, MAF_CHR), (start_col, MAF_START), (end_col, MAF_END),
        (ref_col, MAF_REF),   (alt_col,   MAF_ALT),
    ]:
        if user_col != canonical and user_col in df.columns:
            col_remap[user_col] = canonical
    if col_remap:
        df = df.rename(columns=col_remap)

    missing_cols = [c for c in REQUIRED_MAF_COLS if c not in df.columns]
    if missing_cols:
        rprint(
            f"[bold red]ERROR:[/bold red] Missing required columns: {', '.join(missing_cols)}\n"
            "  Tip: run [cyan]inspect[/cyan] to see column names, "
            "then use --chrom-col / --start-col etc. to remap."
        )
        raise typer.Exit(1)

    for col, default in {gene_col: "UNKNOWN", prot_col: "", cdna_col: "", sample_col: "SAMPLE"}.items():
        if col not in df.columns:
            df[col] = default

    n_total = len(df)
    n_samples_total = df[sample_col].nunique() if sample_col in df.columns else "?"
    console.print(f"  {n_total:,} variants · {n_samples_total} samples (before filtering)")

    # Apply sample filter
    n_filtered_out = 0
    if sample_filter is not None:
        if sample_col not in df.columns:
            rprint(
                f"[bold red]ERROR:[/bold red] Sample column '{sample_col}' not found in MAF. "
                "Cannot apply --samples filter."
            )
            raise typer.Exit(1)

        # Warn about IDs in the filter that don't appear in the MAF at all
        maf_sample_ids = set(df[sample_col].dropna().astype(str).unique())
        unknown_ids = sample_filter - maf_sample_ids
        if unknown_ids:
            console.print(
                f"  [yellow]⚠ {len(unknown_ids)} requested sample ID(s) not found in MAF:[/yellow] "
                + ", ".join(sorted(unknown_ids)[:6])
                + (" …" if len(unknown_ids) > 6 else "")
            )

        before = len(df)
        df = df[df[sample_col].astype(str).isin(sample_filter)].copy()
        n_filtered_out = before - len(df)
        console.print(
            f"  [green]→ {len(df):,} variants kept[/green] · "
            f"{n_filtered_out:,} excluded by sample filter"
        )

        if df.empty:
            rprint("[bold red]ERROR:[/bold red] No variants remain after sample filtering.")
            raise typer.Exit(1)

    # ------------------------------------------------------------------
    # Open reference FASTA and detect chr notation
    # ------------------------------------------------------------------
    console.print(f"[bold]Reference FASTA:[/bold] {ref_genome}  [dim]({genome_build})[/dim]")
    try:
        fa = pysam.FastaFile(str(ref_genome))
    except Exception as exc:
        rprint(f"[bold red]ERROR:[/bold red] Could not open FASTA: {exc}")
        raise typer.Exit(1)

    fasta_has_chr = _detect_fasta_chr_style(fa)
    maf_chr_style = _detect_maf_chr_style(df, MAF_CHR)
    maf_chr_desc  = (
        "chr-prefixed (e.g. chr7)" if maf_chr_style is True
        else "bare (e.g. 7)"        if maf_chr_style is False
        else "unknown / mixed"
    )
    console.print(
        f"  [dim]FASTA contigs: {'chr-prefixed' if fasta_has_chr else 'bare'}  |  "
        f"MAF Chromosome column: {maf_chr_desc}[/dim]"
    )
    if maf_chr_style is not None and maf_chr_style != fasta_has_chr:
        console.print(
            "  [dim]ℹ  Chromosome notation differs between MAF and FASTA — "
            "will be normalised automatically per row.[/dim]"
        )

    # ------------------------------------------------------------------
    # Lazy per-chromosome tabix cache
    # Also tracks chr notation per VCF file (detected on first open).
    # ------------------------------------------------------------------
    tabix_cache:         Dict[str, Optional[object]] = {}
    tabix_has_chr_cache: Dict[str, bool]             = {}

    def _get_tabix(bare: str):
        if bare in tabix_cache:
            return tabix_cache[bare]
        if pop_vcf_dir is None:
            tabix_cache[bare] = None
            return None
        vcf_path = _resolve_vcf_for_chrom(pop_vcf_dir, bare, genome_build)
        if vcf_path is None:
            console.print(
                f"  [yellow]⚠ No VCF for chr{bare} in {pop_vcf_dir.name} — "
                f"masking skipped for this chromosome[/yellow]"
            )
            tabix_cache[bare] = None
            return None
        tbi = Path(str(vcf_path) + ".tbi")
        if not tbi.exists():
            console.print(
                f"  [yellow]⚠ TBI index missing: {vcf_path.name} — "
                f"masking skipped for chr{bare}[/yellow]\n"
                f"    Fix: [cyan]tabix -p vcf {vcf_path}[/cyan]"
            )
            tabix_cache[bare] = None
            return None
        try:
            tbx = pysam.TabixFile(str(vcf_path))
            tabix_cache[bare] = tbx
            tabix_has_chr_cache[bare] = _detect_tabix_chr_style(tbx)
            console.print(
                f"  [dim]Opened {vcf_path.name}  "
                f"(contigs: {'chr-prefixed' if tabix_has_chr_cache[bare] else 'bare'})[/dim]"
            )
            return tbx
        except Exception as exc:
            console.print(f"  [yellow]⚠ Could not open {vcf_path.name}: {exc}[/yellow]")
            tabix_cache[bare] = None
            return None

    if pop_vcf_dir is not None:
        console.print(
            f"[bold]Population VCF dir:[/bold] {pop_vcf_dir}  "
            f"[dim](AF ≥ {af_threshold}, per-chromosome, auto-resolved)[/dim]"
        )
    else:
        console.print("[yellow]⚠ No --pop-vcf-dir — SNP masking skipped.[/yellow]")

    console.print(f"[bold]Flank:[/bold] ±{flank} bp\n")

    # ------------------------------------------------------------------
    # Process variants
    # ------------------------------------------------------------------
    records:        List[str]  = []
    skipped:        int        = 0
    n_masked_total: int        = 0
    summary_rows:   List[dict] = []
    skip_reasons:   List[str]  = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("({task.completed}/{task.total})"),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Processing variants…", total=len(df))

        for row_idx, row in df.iterrows():
            progress.advance(task)

            ref    = str(row[MAF_REF]) if pd.notna(row[MAF_REF]) else "-"
            alt    = str(row[MAF_ALT]) if pd.notna(row[MAF_ALT]) else "-"
            gene   = str(row.get(gene_col,   "UNKNOWN"))
            pann   = str(row.get(prot_col,   ""))
            cann   = str(row.get(cdna_col,   ""))
            sample = str(row.get(sample_col, "SAMPLE"))

            # Normalise chromosome — handles NaN, any case of 'chr' prefix,
            # numeric encodings (23→X, 24→Y, 25/26→MT), mito aliases (M→MT),
            # and validates against VALID_CHROMS.
            bare, chrom_err = _normalise_maf_chrom(row[MAF_CHR])
            if chrom_err:
                skip_reasons.append(f"row {row_idx} {gene} — {chrom_err}")
                skipped += 1
                continue
            # Preserve the original prefix style for display/header only
            raw_chrom_str = str(row[MAF_CHR]).strip()
            chrom = raw_chrom_str if raw_chrom_str.upper().startswith("CHR") else bare

            try:
                start = int(float(row[MAF_START]))
                end   = int(float(row[MAF_END]))
            except (ValueError, TypeError):
                skip_reasons.append(
                    f"row {row_idx}: non-numeric position "
                    f"(start={row[MAF_START]!r}, end={row[MAF_END]!r})"
                )
                skipped += 1
                continue

            coord_err = _validate_coordinates(start, end)
            if coord_err:
                reason = f"row {row_idx} {gene} {chrom}:{start} — {coord_err}"
                if skip_invalid:
                    skip_reasons.append(reason)
                    skipped += 1
                    continue
                rprint(f"[bold red]ERROR:[/bold red] {reason}")
                raise typer.Exit(1)

            if not _validate_allele(ref) or not _validate_allele(alt):
                reason = (
                    f"row {row_idx} {gene} {chrom}:{start} — "
                    f"invalid allele (ref={ref!r}, alt={alt!r})"
                )
                if skip_invalid:
                    skip_reasons.append(reason)
                    skipped += 1
                    continue
                rprint(f"[bold red]ERROR:[/bold red] {reason}")
                raise typer.Exit(1)

            # Ensure tabix is open (lazy, cached); notation detected on open
            tabix = _get_tabix(bare)

            try:
                fl, fr, ml, mr = _fetch_flanks(
                    fa, tabix_cache,
                    fasta_has_chr, tabix_has_chr_cache,
                    bare, start, end,
                    flank, af_threshold,
                )
            except Exception as exc:
                skip_reasons.append(f"row {row_idx} {gene} {chrom}:{start} — fetch error: {exc}")
                skipped += 1
                continue

            if uppercase:
                fl, fr, ml, mr = fl.upper(), fr.upper(), ml.upper(), mr.upper()
                ref = ref.upper()
                alt = alt.upper()

            n_masked = ml.count("N") + mr.count("N")
            n_masked_total += n_masked

            header_base = f"{_safe(sample)}__{_safe(gene)}__{_safe(pann)}__{_safe(cann)}"
            records.append(f">{header_base}\n{fl}[{ref}/{alt}]{fr}\n")
            records.append(f">Masked__{header_base}\n{ml}[{ref}/{alt}]{mr}\n")

            summary_rows.append({
                "Sample":   sample[:32],
                "Gene":     gene,
                "Chrom":    chrom,
                "Pos":      f"{start}-{end}",
                "Ref":      ref[:8],
                "Alt":      alt[:8],
                "Left bp":  len(fl),
                "Right bp": len(fr),
                "N masked": n_masked,
            })

    # ------------------------------------------------------------------
    # Close handles
    # ------------------------------------------------------------------
    fa.close()
    for tbx in tabix_cache.values():
        if tbx is not None:
            tbx.close()

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    try:
        output.write_text("".join(records))
    except OSError as exc:
        rprint(f"[bold red]ERROR:[/bold red] Could not write output: {exc}")
        raise typer.Exit(1)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    console.rule("[bold green]Results[/bold green]")

    if summary_rows:
        table = Table(show_header=True, header_style="bold cyan", show_lines=False, highlight=True)
        for col in ("Sample", "Gene", "Chrom", "Pos", "Ref", "Alt", "Left bp", "Right bp", "N masked"):
            table.add_column(col, no_wrap=(col in ("Sample", "Gene", "Chrom", "Pos")))
        for r in summary_rows:
            n = r["N masked"]
            n_str = f"[yellow]{n}[/yellow]" if n > 0 else "[dim]0[/dim]"
            table.add_row(
                r["Sample"], r["Gene"], r["Chrom"], r["Pos"],
                r["Ref"], r["Alt"],
                str(r["Left bp"]), str(r["Right bp"]), n_str,
            )
        console.print(table)

    if skip_reasons:
        console.print(f"\n[bold yellow]Skipped {skipped} variants:[/bold yellow]")
        for reason in skip_reasons[:20]:
            console.print(f"  • {reason}")
        if len(skip_reasons) > 20:
            console.print(f"  [dim]… and {len(skip_reasons) - 20} more[/dim]")

    filter_line = (
        f"[bold]Excluded by filter:[/bold] {n_filtered_out:,}\n"
        if sample_filter is not None else ""
    )
    console.print(
        f"\n[bold]Total in MAF:[/bold]     {n_total:>6,}\n"
        + filter_line +
        f"[bold]Processed:[/bold]        {len(summary_rows):>6,}\n"
        f"[bold]Skipped (errors):[/bold] {skipped:>6,}\n"
        f"[bold]Bases masked:[/bold]     {n_masked_total:>6,}\n"
        f"[bold]FASTA records:[/bold]    {len(records):>6,} [dim](2 per variant: raw + masked)[/dim]\n"
        f"[bold]Output:[/bold] [cyan]{output.resolve()}[/cyan]\n"
        f"[dim]Elapsed: {elapsed:.1f}s[/dim]"
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app()
