"""`vflank fusion` — structural-variant junction extraction (simple TSV input).

Reads a breakpoint table (chr1 pos1 str1 chr2 pos2 str2; columns matched by
name) and writes one FASTA record per fusion: the chimeric junction sequence a
ddPCR probe spans. VCF/BND input is a later phase.
"""

from __future__ import annotations

import time
from pathlib import Path

import typer
from rich.table import Table

from ..core.fusion import build_junction
from ..errors import VflankError
from ..io import breakpoints as bp_io
from ..io import fasta as fasta_io
from ..io.breakpoints import SvColumns
from ..io.reference import ReferenceFasta
from ..logging import console

app = typer.Typer(no_args_is_help=True)


@app.command()
def run(
    sv_file: Path = typer.Argument(
        ..., exists=True, help="Breakpoint TSV (chr1 pos1 str1 chr2 pos2 str2)."
    ),
    ref_genome: Path = typer.Option(
        ..., "--ref-genome", "-r", help="Indexed reference FASTA (.fai required)."
    ),
    genome_build: str = typer.Option("hg19", "--genome-build", "-g", help="hg19 or hg38."),
    flank: int = typer.Option(
        200, "--flank", "-f", min=1, max=10_000,
        help="Bases taken from each partner (junction is up to 2x this).",
    ),
    output: Path = typer.Option(
        Path("fusion_junctions.fasta"), "--output", "-o", help="Output FASTA file."
    ),
    chr1_col: str = typer.Option(SvColumns.chr1, "--chr1-col"),
    pos1_col: str = typer.Option(SvColumns.pos1, "--pos1-col"),
    str1_col: str = typer.Option(SvColumns.str1, "--str1-col"),
    chr2_col: str = typer.Option(SvColumns.chr2, "--chr2-col"),
    pos2_col: str = typer.Option(SvColumns.pos2, "--pos2-col"),
    str2_col: str = typer.Option(SvColumns.str2, "--str2-col"),
    name_col: str = typer.Option(SvColumns.name, "--name-col"),
    sample_col: str = typer.Option(SvColumns.sample, "--sample-col"),
):
    """Build fusion-junction sequences for a breakpoint table and write a FASTA."""
    cols = SvColumns(
        chr1_col, pos1_col, str1_col, chr2_col, pos2_col, str2_col, name_col, sample_col
    )
    try:
        _run(sv_file, ref_genome, genome_build, flank, output, cols)
    except VflankError as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        raise typer.Exit(1) from exc


def _run(sv_file, ref_genome, genome_build, flank, output, cols: SvColumns):
    t0 = time.time()
    console.rule("[bold blue]vflank fusion run[/bold blue]")
    if genome_build not in ("hg19", "hg38"):
        raise VflankError(f"--genome-build must be 'hg19' or 'hg38', got '{genome_build}'")

    console.print(f"[bold]Loading breakpoints:[/bold] {sv_file}")
    df = bp_io.load_sv_table(sv_file, cols)
    console.print(f"  {len(df):,} fusion(s)")

    reference = ReferenceFasta(ref_genome)
    console.print(f"[bold]Reference:[/bold] {ref_genome}  [dim]({genome_build})[/dim]")
    build_warn = reference.check_build(genome_build)
    if build_warn:
        console.print(f"  [bold yellow]⚠ {build_warn}[/bold yellow]")
    console.print(f"[bold]Flank:[/bold] {flank} bp/partner (junction ≤ {2 * flank} bp)\n")

    records: list[str] = []
    skipped = 0
    skip_reasons: list[str] = []
    summary_rows: list[dict] = []

    for row_idx, row in df.iterrows():
        fusion, reason = bp_io.parse_fusion_row(row, cols)
        if reason is not None:
            skip_reasons.append(f"row {row_idx} — {reason}")
            skipped += 1
            continue
        try:
            jr = build_junction(reference, fusion, flank)
        except Exception as exc:  # noqa: BLE001
            skip_reasons.append(f"row {row_idx} {fusion.name} — junction error: {exc}")
            skipped += 1
            continue

        label = fasta_io.safe_header(fusion.name or "fusion")
        bp = f"{fusion.bp1.chrom}_{fusion.bp1.pos}_{fusion.bp1.strand}__" \
             f"{fusion.bp2.chrom}_{fusion.bp2.pos}_{fusion.bp2.strand}"
        prefix = f"{fasta_io.safe_header(fusion.sample)}__" if fusion.sample else ""
        header = f"{prefix}{label}__{bp}__j{jr.junction_index}"
        records.append(f">{header}\n{jr.sequence}\n")

        truncated = len(jr.sequence) < 2 * flank
        summary_rows.append({
            "Name": fusion.name or ".", "BP1": f"{fusion.bp1.chrom}:{fusion.bp1.pos}",
            "BP2": f"{fusion.bp2.chrom}:{fusion.bp2.pos}",
            "Len": len(jr.sequence), "Junction": jr.junction_index,
            "Trunc": truncated,
        })

    reference.close()
    fasta_io.write_fasta(output, records)

    console.rule("[bold green]Results[/bold green]")
    if summary_rows:
        table = Table(show_header=True, header_style="bold cyan")
        for col in ("Name", "BP1", "BP2", "Len", "Junction", "Trunc"):
            table.add_column(col)
        for r in summary_rows[:50]:
            table.add_row(
                r["Name"], r["BP1"], r["BP2"], str(r["Len"]),
                str(r["Junction"]), "[yellow]yes[/yellow]" if r["Trunc"] else "no",
            )
        console.print(table)

    if skip_reasons:
        console.print(f"\n[bold yellow]Skipped {skipped}:[/bold yellow]")
        for reason in skip_reasons[:20]:
            console.print(f"  • {reason}")

    console.print(
        f"\n[bold]Fusions:[/bold]   {len(df):>6,}\n"
        f"[bold]Junctions:[/bold] {len(records):>6,}\n"
        f"[bold]Skipped:[/bold]   {skipped:>6,}\n"
        f"[bold]Output:[/bold] [cyan]{output.resolve()}[/cyan]  "
        f"[dim]({time.time() - t0:.1f}s)[/dim]"
    )
