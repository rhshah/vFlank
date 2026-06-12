"""`vflank small` — small-variant (SNP/indel) flank extraction and masking."""

from __future__ import annotations

import time
from collections import Counter
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

from ..core.chrom import detect_series_chr_style, normalise_chrom
from ..core.flanks import ReferenceFlankSource
from ..core.popfreq import GnomadStore, build_chrom_vcf_map, kinds_for
from ..core.popfreq_api import GnomadApiSource
from ..core.skips import categorize_skip
from ..errors import VflankError
from ..io import fasta as fasta_io
from ..io import maf as maf_io
from ..io import report as report_io
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
        "hg19", "--genome-build", "-g",
        help="hg19 (GRCh37, gnomAD v2.1.1) or hg38 (GRCh38, gnomAD v4).",
    ),
    flank: int = typer.Option(
        200, "--flank", "-f", min=1, max=10_000, help="Bases on each side of the variant."
    ),
    af_threshold: float = typer.Option(
        0.001, "--af-threshold", min=0.0, max=1.0, help="Min population AF to mask a SNP."
    ),
    pop_data: str = typer.Option(
        "genome", "--pop-data",
        help="gnomAD data to mask against: genome (default), exome, or both (union).",
    ),
    pop_source: str = typer.Option(
        "vcf", "--pop-source",
        help="Masking backend: vcf (local gnomAD VCFs) or api (gnomAD GraphQL, no download).",
    ),
    output: Path = typer.Option(
        Path("flanking_sequences.fasta"), "--output", "-o", help="Output FASTA file."
    ),
    report: Path | None = typer.Option(
        None, "--report",
        help="Write a per-variant TSV run report (stats + table) to this path.",
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
    dedup: bool = typer.Option(
        True, "--dedup/--no-dedup",
        help="Emit one record per unique variant (CHR_POS_REF_ALT), collapsing samples.",
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
            pop_data, pop_source, output, report, samples, samples_file,
            MafColumns(chrom_col, start_col, end_col, ref_col, alt_col,
                       gene_col, prot_col, cdna_col, sample_col),
            uppercase, dedup,
        )
    except VflankError as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        raise typer.Exit(1) from exc


def _run(maf_file, ref_genome, pop_vcf_dir, genome_build, flank, af_threshold,
         pop_data, pop_source, output, report, samples, samples_file,
         cols: MafColumns, uppercase: bool, dedup: bool):
    t0 = time.time()
    console.rule("[bold blue]vflank small run[/bold blue]")

    if genome_build not in ("hg19", "hg38"):
        raise VflankError(f"--genome-build must be 'hg19' or 'hg38', got '{genome_build}'")
    if pop_data not in ("genome", "exome", "both"):
        raise VflankError(f"--pop-data must be 'genome', 'exome', or 'both', got '{pop_data}'")
    if pop_source not in ("vcf", "api"):
        raise VflankError(f"--pop-source must be 'vcf' or 'api', got '{pop_source}'")
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

    # --- Masking source (optional) ---
    gnomad = None
    if pop_source == "api":
        gnomad = GnomadApiSource(genome_build, pop_data)
        if pop_vcf_dir is not None:
            console.print("  [yellow]⚠ --pop-vcf-dir ignored with --pop-source api[/yellow]")
        if len(df) > 50:
            console.print(
                f"  [yellow]⚠ {len(df):,} variants over the rate-limited API "
                f"(~10 req/min) — consider --pop-source vcf for bulk.[/yellow]"
            )
        console.print(
            f"[bold]Masking:[/bold] gnomAD API  "
            f"[dim](pop-data={pop_data}, dataset={gnomad.dataset}, AF ≥ {af_threshold})[/dim]"
        )
    elif pop_vcf_dir is not None:
        gnomad = GnomadStore(pop_vcf_dir, genome_build, pop_data)
        # Fail fast if a requested data kind is wholly absent for the MAF's
        # chromosomes (no silent fall-back to genome-only when exome/both asked).
        maf_chroms = {
            b for b, _err in (normalise_chrom(c) for c in df[MAF_CHR].dropna().unique()) if b
        }
        if maf_chroms:
            gnomad.preflight(sorted(maf_chroms))
        console.print(
            f"[bold]Population VCF dir:[/bold] {pop_vcf_dir}  "
            f"[dim](pop-data={pop_data}, AF ≥ {af_threshold})[/dim]"
        )
    else:
        console.print(
            "[yellow]⚠ No masking source (--pop-source vcf without --pop-vcf-dir).[/yellow]"
        )
    console.print(f"[bold]Flank:[/bold] ±{flank} bp\n")

    source = ReferenceFlankSource(reference, gnomad, flank=flank, af_threshold=af_threshold)

    # --- Process variants ---
    records: list[str] = []
    skipped = 0
    n_duplicate = 0
    n_masked_total = 0
    seen_variants: set[tuple] = set()
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

            # Collapse the same variant seen across multiple samples: the flank
            # and mask are sample-independent here, so one record per unique
            # CHR_POS_REF_ALT. (Consensus modes C/D will key on variant+sample.)
            if dedup:
                key = (variant.chrom, variant.start, variant.end,
                       variant.ref.upper(), variant.alt.upper())
                if key in seen_variants:
                    n_duplicate += 1
                    continue
                seen_variants.add(key)

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
            truncated = len(fr.left) < exp_left or len(fr.right) < flank
            if truncated:
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
                "Chrom": variant.raw_chrom, "Start": variant.start, "End": variant.end,
                "Ref": ref[:8], "Alt": alt[:8],
                "LeftLen": len(fr.left), "RightLen": len(fr.right),
                "NMasked": fr.n_masked, "Truncated": truncated,
            })

    reference.close()
    api_requests = getattr(gnomad, "request_count", None) if gnomad is not None else None
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
            n_str = f"[yellow]{r['NMasked']}[/yellow]" if r["NMasked"] else "[dim]0[/dim]"
            table.add_row(
                r["Sample"], r["Gene"], r["Chrom"], f"{r['Start']}-{r['End']}",
                r["Ref"], r["Alt"], str(r["LeftLen"]), str(r["RightLen"]), n_str,
            )
        console.print(table)
        if len(summary_rows) > 50:
            console.print(f"  [dim]… and {len(summary_rows) - 50} more[/dim]")

    # Categorise skips so a large, uniform skip set (e.g. 91 missing chromosomes)
    # reads as one line per category instead of a wall of identical messages.
    skip_breakdown = Counter(categorize_skip(r) for r in skip_reasons)
    if skip_reasons:
        console.print(f"\n[bold yellow]Skipped {skipped} — by reason:[/bold yellow]")
        for category, count in skip_breakdown.most_common():
            console.print(f"  [yellow]{count:>5,}[/yellow]  {category}")
        console.print("  [dim]examples:[/dim]")
        for reason in skip_reasons[:5]:
            console.print(f"    • {reason}")
        if len(skip_reasons) > 5:
            console.print(f"    [dim]… {len(skip_reasons) - 5} more (see --report)[/dim]")

    if flank_warnings:
        console.print(
            f"\n[bold yellow]Truncated flanks {len(flank_warnings)} "
            "(emitted, but shorter than requested):[/bold yellow]"
        )
        for reason in flank_warnings[:10]:
            console.print(f"  • {reason}")
        if len(flank_warnings) > 10:
            console.print(f"  [dim]… and {len(flank_warnings) - 10} more[/dim]")

    n_truncated = len(flank_warnings)
    truncated_line = (
        f"[bold]Truncated flanks:[/bold] {n_truncated:>6,}\n" if n_truncated else ""
    )
    api_line = (
        f"[bold]API requests:[/bold]     {api_requests:>6,}\n" if api_requests is not None else ""
    )
    dup_line = (
        f"[bold]Dup. collapsed:[/bold]   {n_duplicate:>6,} [dim](other samples)[/dim]\n"
        if n_duplicate else ""
    )
    console.print(
        f"\n[bold]Total in MAF:[/bold]     {n_total:>6,}\n"
        f"[bold]Processed:[/bold]        {len(summary_rows):>6,}\n"
        f"[bold]Skipped:[/bold]          {skipped:>6,}\n"
        + dup_line + truncated_line +
        f"[bold]Bases masked:[/bold]     {n_masked_total:>6,}\n"
        + api_line +
        f"[bold]FASTA records:[/bold]    {len(records):>6,} [dim](2 per variant)[/dim]\n"
        f"[bold]Output:[/bold] [cyan]{output.resolve()}[/cyan]  [dim]({elapsed:.1f}s)[/dim]"
    )

    # --- Optional machine-readable report ---
    if report is not None:
        stats = {
            "total_in_maf": n_total,
            "processed": len(summary_rows),
            "skipped": skipped,
            "duplicates_collapsed": n_duplicate,
            "truncated_flanks": n_truncated,
            "bases_masked": n_masked_total,
            "fasta_records": len(records),
            "pop_source": pop_source,
            "pop_data": pop_data,
        }
        if api_requests is not None:
            stats["api_requests"] = api_requests
        report_io.write_report(report, summary_rows, stats, dict(skip_breakdown))
        console.print(f"[bold]Report:[/bold] [cyan]{report.resolve()}[/cyan]")


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
    genome_build: str = typer.Option("hg19", "--genome-build", "-g", help="hg19 or hg38"),
    pop_data: str = typer.Option(
        "genome", "--pop-data", help="Coverage to check: genome, exome, or both."
    ),
):
    """Show which per-chromosome VCFs were found (and which are missing) per data kind."""
    if pop_data not in ("genome", "exome", "both"):
        console.print(
            f"[bold red]ERROR:[/bold red] --pop-data must be genome|exome|both, got '{pop_data}'"
        )
        raise typer.Exit(1)

    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    kinds = kinds_for(pop_data)
    kind_maps = {k: build_chrom_vcf_map(pop_vcf_dir, genome_build, chroms, k) for k in kinds}

    table = Table(show_header=True, header_style="bold cyan",
                  title=f"VCFs in {pop_vcf_dir}  [dim]({genome_build}, {pop_data})[/dim]")
    table.add_column("Chrom", width=6)
    for kind in kinds:
        table.add_column(kind, overflow="fold")
        table.add_column(f"{kind} TBI", justify="center", width=10)

    found = {k: 0 for k in kinds}
    for chrom in chroms:
        cells: list[str] = [chrom]
        for kind in kinds:
            path = kind_maps[kind].get(chrom)
            if path:
                tbi_ok = Path(str(path) + ".tbi").exists()
                cells.append(path.name)
                cells.append("[green]✓[/green]" if tbi_ok else "[red]✗[/red]")
                found[kind] += 1
            else:
                cells.append("[dim]not found[/dim]")
                cells.append("[dim]—[/dim]")
        table.add_row(*cells)
    console.print(table)
    for kind in kinds:
        console.print(f"{found[kind]}/{len(chroms)} chromosomes have a matching {kind} VCF.")
