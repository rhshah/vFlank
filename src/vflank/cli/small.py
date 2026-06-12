"""`vflank small` — small-variant (SNP/indel) flank extraction and masking."""

from __future__ import annotations

import time
from pathlib import Path

import typer
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from ..core.chrom import detect_series_chr_style
from ..core.flanks import ReferenceFlankSource
from ..core.popfreq import GnomadStore, build_chrom_vcf_map
from ..errors import VflankError
from ..io import fasta as fasta_io
from ..io import maf as maf_io
from ..io.maf import MAF_CHR, MAF_SAMPLE, REQUIRED_MAF_COLS, MafColumns
from ..io.reference import ReferenceFasta
from ..logging import console

app = typer.Typer(no_args_is_help=True)


def _load_sample_filter(
    samples: str | None, samples_file: Path | None
) -> set[str] | None:
    ids: set[str] = set()
    if samples:
        ids |= {s.strip() for s in samples.split(",") if s.strip()}
    if samples_file is not None:
        if not samples_file.exists():
            raise VflankError(f"--samples-file not found: {samples_file}")
        for line in samples_file.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                ids.add(line)
    if not ids and (samples is not None or samples_file is not None):
        raise VflankError("--samples / --samples-file produced an empty ID set.")
    return ids or None


@app.command()
def run(
    maf_file: Path = typer.Argument(
        ..., help="Input MAF (tab-separated, TCGA/MSK).", exists=True
    ),
    ref_genome: Path = typer.Option(
        ..., "--ref-genome", "-r", help="Indexed reference FASTA (.fai required)."
    ),
    pop_vcf_dir: Path | None = typer.Option(
        None, "--pop-vcf-dir", "-d",
        help="Directory of per-chromosome gnomAD VCF bgz files. Omit to skip masking.",
    ),
    genome_build: str = typer.Option(
        "hg38", "--genome-build", "-g", help="hg38 (GRCh38) or hg19 (GRCh37)."
    ),
    flank: int = typer.Option(
        200, "--flank", "-f", min=1, max=10_000, help="Bases on each side of the variant."
    ),
    af_threshold: float = typer.Option(
        0.001, "--af-threshold", min=0.0, max=1.0, help="Min population AF to mask a SNP."
    ),
    output: Path = typer.Option(
        Path("flanking_sequences.fasta"), "--output", "-o", help="Output FASTA file."
    ),
    samples: str | None = typer.Option(
        None, "--samples", "-s",
        help="Comma-separated Tumor_Sample_Barcode IDs to include.",
    ),
    samples_file: Path | None = typer.Option(
        None, "--samples-file",
        help="File of sample IDs, one per line (# comments allowed).",
    ),
    chrom_col: str = typer.Option(MafColumns.chrom, "--chrom-col"),
    start_col: str = typer.Option(MafColumns.start, "--start-col"),
    end_col: str = typer.Option(MafColumns.end, "--end-col"),
    ref_col: str = typer.Option(MafColumns.ref, "--ref-col"),
    alt_col: str = typer.Option(MafColumns.alt, "--alt-col"),
    gene_col: str = typer.Option(MafColumns.gene, "--gene-col"),
    prot_col: str = typer.Option(MafColumns.protein, "--prot-col"),
    cdna_col: str = typer.Option(MafColumns.cdna, "--cdna-col"),
    sample_col: str = typer.Option(MafColumns.sample, "--sample-col"),
    uppercase: bool = typer.Option(
        True, "--uppercase/--no-uppercase", help="Uppercase flanking sequences."
    ),
):
    """Extract flanking sequences for every variant in a MAF and write a FASTA.

    Each variant yields a raw record and a Masked record (common population SNPs
    with AF >= --af-threshold replaced by 'N'). Chromosome notation is
    auto-detected from the FASTA and VCFs.
    """
    try:
        _run(
            maf_file, ref_genome, pop_vcf_dir, genome_build, flank, af_threshold,
            output, samples, samples_file,
            MafColumns(chrom_col, start_col, end_col, ref_col, alt_col,
                       gene_col, prot_col, cdna_col, sample_col),
            uppercase,
        )
    except VflankError as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        raise typer.Exit(1) from exc


def _run(maf_file, ref_genome, pop_vcf_dir, genome_build, flank, af_threshold,
         output, samples, samples_file, cols: MafColumns, uppercase: bool):
    t0 = time.time()
    console.rule("[bold blue]vflank small run[/bold blue]")

    if genome_build not in ("hg19", "hg38"):
        raise VflankError(f"--genome-build must be 'hg19' or 'hg38', got '{genome_build}'")
    if pop_vcf_dir is not None and not pop_vcf_dir.is_dir():
        raise VflankError(f"--pop-vcf-dir is not a directory: {pop_vcf_dir}")

    sample_filter = _load_sample_filter(samples, samples_file)
    if sample_filter is not None:
        console.print(f"[bold]Sample filter:[/bold] {len(sample_filter)} ID(s)")

    # --- Load + filter MAF ---
    console.print(f"[bold]Loading MAF:[/bold] {maf_file}")
    df = maf_io.load_maf(maf_file, cols)
    n_total = len(df)
    console.print(f"  {n_total:,} variants before filtering")

    n_filtered_out = 0
    if sample_filter is not None:
        if cols.sample not in df.columns:
            raise VflankError(f"Sample column '{cols.sample}' not found; cannot filter.")
        maf_ids = set(df[cols.sample].dropna().astype(str).unique())
        unknown = sample_filter - maf_ids
        if unknown:
            console.print(f"  [yellow]⚠ {len(unknown)} requested ID(s) not in MAF[/yellow]")
        before = len(df)
        df = df[df[cols.sample].astype(str).isin(sample_filter)].copy()
        n_filtered_out = before - len(df)
        console.print(f"  [green]→ {len(df):,} kept[/green] · {n_filtered_out:,} excluded")
        if df.empty:
            raise VflankError("No variants remain after sample filtering.")

    # --- Reference FASTA + build guard ---
    reference = ReferenceFasta(ref_genome)
    console.print(f"[bold]Reference:[/bold] {ref_genome}  [dim]({genome_build})[/dim]")
    build_warn = reference.check_build(genome_build)
    if build_warn:
        console.print(f"  [bold yellow]⚠ {build_warn}[/bold yellow]")

    maf_style = detect_series_chr_style(df[MAF_CHR].head(20)) if MAF_CHR in df.columns else None
    console.print(
        f"  [dim]FASTA contigs: {'chr-prefixed' if reference.has_chr else 'bare'} · "
        f"MAF: {'chr-prefixed' if maf_style else 'bare' if maf_style is False else 'unknown'}[/dim]"
    )

    # --- Population VCF store (optional) ---
    gnomad = None
    if pop_vcf_dir is not None:
        gnomad = GnomadStore(pop_vcf_dir, genome_build)
        console.print(
            f"[bold]Population VCF dir:[/bold] {pop_vcf_dir}  [dim](AF ≥ {af_threshold})[/dim]"
        )
    else:
        console.print("[yellow]⚠ No --pop-vcf-dir — SNP masking skipped.[/yellow]")
    console.print(f"[bold]Flank:[/bold] ±{flank} bp\n")

    source = ReferenceFlankSource(reference, gnomad, flank=flank, af_threshold=af_threshold)

    # --- Process variants ---
    records: list[str] = []
    skipped = 0
    n_masked_total = 0
    summary_rows: list[dict] = []
    skip_reasons: list[str] = []
    flank_warnings: list[str] = []

    with Progress(
        SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
        BarColumn(), TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("({task.completed}/{task.total})"), TimeElapsedColumn(),
        console=console, transient=True,
    ) as progress:
        task = progress.add_task("Processing variants…", total=len(df))
        for row_idx, row in df.iterrows():
            progress.advance(task)

            variant, reason = maf_io.parse_variant_row(row, cols)
            if reason is not None:
                skip_reasons.append(f"row {row_idx} {reason}")
                skipped += 1
                continue

            try:
                fr = source.fetch(variant)
            except Exception as exc:  # noqa: BLE001
                skip_reasons.append(
                    f"row {row_idx} {variant.gene} {variant.raw_chrom}:{variant.start} "
                    f"— fetch error: {exc}"
                )
                skipped += 1
                continue

            # Flag silently-truncated flanks. pysam returns a short string when a
            # window runs off a contig end rather than raising. A short left flank
            # near position 1 is expected; anything shorter than requested on the
            # right (or shorter than the position allows on the left) is a contig
            # boundary the user should know about — the record is still emitted.
            exp_left = min(flank, variant.start - 1)
            if len(fr.left) < exp_left or len(fr.right) < flank:
                flank_warnings.append(
                    f"row {row_idx} {variant.gene} {variant.raw_chrom}:{variant.start} — "
                    f"flank truncated near contig boundary "
                    f"(left {len(fr.left)}/{exp_left}, right {len(fr.right)}/{flank})"
                )

            ref, alt = variant.ref, variant.alt
            if uppercase:
                fr = fr.upper()
                ref, alt = ref.upper(), alt.upper()

            n_masked_total += fr.n_masked
            records.extend(fasta_io.format_records(variant, fr, ref, alt))
            summary_rows.append({
                "Sample": variant.sample[:32], "Gene": variant.gene,
                "Chrom": variant.raw_chrom, "Pos": f"{variant.start}-{variant.end}",
                "Ref": ref[:8], "Alt": alt[:8],
                "Left": len(fr.left), "Right": len(fr.right), "N": fr.n_masked,
            })

    reference.close()
    if gnomad is not None:
        gnomad.close()

    fasta_io.write_fasta(output, records)

    # --- Summary ---
    elapsed = time.time() - t0
    console.rule("[bold green]Results[/bold green]")
    if summary_rows:
        table = Table(show_header=True, header_style="bold cyan", highlight=True)
        for col in ("Sample", "Gene", "Chrom", "Pos", "Ref", "Alt", "Left", "Right", "N"):
            table.add_column(col, no_wrap=col in ("Sample", "Gene", "Chrom", "Pos"))
        for r in summary_rows[:50]:
            n_str = f"[yellow]{r['N']}[/yellow]" if r["N"] else "[dim]0[/dim]"
            table.add_row(r["Sample"], r["Gene"], r["Chrom"], r["Pos"], r["Ref"],
                          r["Alt"], str(r["Left"]), str(r["Right"]), n_str)
        console.print(table)
        if len(summary_rows) > 50:
            console.print(f"  [dim]… and {len(summary_rows) - 50} more[/dim]")

    if skip_reasons:
        console.print(f"\n[bold yellow]Skipped {skipped}:[/bold yellow]")
        for reason in skip_reasons[:20]:
            console.print(f"  • {reason}")
        if len(skip_reasons) > 20:
            console.print(f"  [dim]… and {len(skip_reasons) - 20} more[/dim]")

    if flank_warnings:
        console.print(
            f"\n[bold yellow]Truncated flanks {len(flank_warnings)} "
            "(emitted, but shorter than requested):[/bold yellow]"
        )
        for reason in flank_warnings[:20]:
            console.print(f"  • {reason}")
        if len(flank_warnings) > 20:
            console.print(f"  [dim]… and {len(flank_warnings) - 20} more[/dim]")

    truncated_line = (
        f"[bold]Truncated flanks:[/bold] {len(flank_warnings):>6,}\n" if flank_warnings else ""
    )
    console.print(
        f"\n[bold]Total in MAF:[/bold]     {n_total:>6,}\n"
        f"[bold]Processed:[/bold]        {len(summary_rows):>6,}\n"
        f"[bold]Skipped:[/bold]          {skipped:>6,}\n"
        + truncated_line +
        f"[bold]Bases masked:[/bold]     {n_masked_total:>6,}\n"
        f"[bold]FASTA records:[/bold]    {len(records):>6,} [dim](2 per variant)[/dim]\n"
        f"[bold]Output:[/bold] [cyan]{output.resolve()}[/cyan]  [dim]({elapsed:.1f}s)[/dim]"
    )


@app.command()
def inspect(
    maf_file: Path = typer.Argument(..., help="MAF file to preview", exists=True),
    n_rows: int = typer.Option(3, "--rows", "-n", min=1, max=20, help="Data rows to show."),
):
    """Print column names and a data preview without running any analysis."""
    df = maf_io.read_maf(maf_file, n_rows=n_rows)
    console.print(Panel(
        f"[bold cyan]{maf_file.name}[/bold cyan]  "
        f"[dim]{len(df.columns)} cols · {len(df)} rows[/dim]",
        title="MAF Inspect",
    ))
    table = Table(show_header=True, header_style="bold magenta", show_lines=True)
    table.add_column("#", style="dim", width=4)
    table.add_column("Column")
    table.add_column("Req?", justify="center", width=6)
    for i in range(min(n_rows, len(df))):
        table.add_column(f"Row {i+1}", overflow="fold")
    for i, col in enumerate(df.columns, 1):
        req = "[green]✓[/green]" if col in REQUIRED_MAF_COLS else ""
        vals = [str(df[col].iloc[j]) for j in range(min(n_rows, len(df)))]
        table.add_row(str(i), col, req, *vals)
    console.print(table)

    missing = [c for c in REQUIRED_MAF_COLS if c not in df.columns]
    if missing:
        console.print(f"\n[bold red]Missing required columns:[/bold red] {', '.join(missing)}")
    else:
        console.print("\n[bold green]✓ All required columns present.[/bold green]")

    if MAF_SAMPLE in df.columns:
        samples = df[MAF_SAMPLE].dropna().unique().tolist()
        console.print(
            f"\n[bold]Samples ({len(samples)}):[/bold] "
            + ", ".join(map(str, samples[:10]))
        )


@app.command("list-vcf")
def list_vcf(
    pop_vcf_dir: Path = typer.Argument(
        ..., help="Directory of gnomAD per-chromosome VCF bgz files.", exists=True
    ),
    genome_build: str = typer.Option("hg38", "--genome-build", "-g", help="hg38 or hg19"),
):
    """Show which per-chromosome VCFs were found (and which are missing)."""
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    vcf_map = build_chrom_vcf_map(pop_vcf_dir, genome_build, chroms)
    table = Table(show_header=True, header_style="bold cyan",
                  title=f"VCFs in {pop_vcf_dir}  [dim]({genome_build})[/dim]")
    table.add_column("Chrom", width=8)
    table.add_column("File found", overflow="fold")
    table.add_column("TBI", justify="center", width=10)
    found = 0
    for chrom in chroms:
        path = vcf_map.get(chrom)
        if path:
            tbi = Path(str(path) + ".tbi")
            tbi_str = "[green]✓[/green]" if tbi.exists() else "[red]✗[/red]"
            table.add_row(chrom, path.name, tbi_str)
            found += 1
        else:
            table.add_row(chrom, "[dim]not found[/dim]", "[dim]—[/dim]")
    console.print(table)
    console.print(f"\n{found}/{len(chroms)} chromosomes have a matching VCF.")
