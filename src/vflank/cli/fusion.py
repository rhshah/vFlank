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

from .. import __version__
from ..core.chrom import normalise_chrom
from ..core.consensus import BamConsensusSource
from ..core.fusion import build_junction
from ..core.popfreq_api import dataset_for_build
from ..errors import ConsensusError, VflankError
from ..io import breakpoints as bp_io
from ..io import fasta as fasta_io
from ..io.breakpoints import SvColumns
from ..io.reference import ReferenceFasta
from ..logging import console, get_logger
from ._bam import build_consensus_policy, load_bam_resolver
from ._masking import make_pop_source, validate_pop_options
from ._ui import echo_parameters

app = typer.Typer(no_args_is_help=True)
log = get_logger()


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
    pop_vcf_dir: Path | None = typer.Option(
        None, "--pop-vcf-dir", "-d",
        help="Directory of gnomAD VCFs to mask junction flanks. Omit to skip masking.",
    ),
    pop_data: str = typer.Option(
        "genome", "--pop-data", help="gnomAD data to mask against: genome, exome, or both."
    ),
    pop_source: str = typer.Option(
        "vcf", "--pop-source", help="Masking backend: vcf or api (no download)."
    ),
    af_threshold: float = typer.Option(
        0.001, "--af-threshold", min=0.0, max=1.0, help="Min population AF to mask a SNP."
    ),
    output: Path = typer.Option(
        Path("fusion_junctions.fasta"), "--output", "-o", help="Output FASTA file."
    ),
    bam: Path | None = typer.Option(
        None, "--bam", help="Single-sample BAM for patient consensus of the junction flanks."
    ),
    bam_map: Path | None = typer.Option(
        None, "--bam-map", help="TSV (sample<TAB>bam_path) for per-fusion consensus."
    ),
    bam_min_depth: int = typer.Option(20, "--bam-min-depth"),
    bam_call_fract: float = typer.Option(0.9, "--bam-call-fract"),
    bam_het_char: str = typer.Option("N", "--bam-het-char", help="Het output: N or iupac."),
    bam_lowcov: str = typer.Option(
        "gnomad", "--bam-lowcov", help="Low-coverage base: n | reference | gnomad."
    ),
    bam_min_baseq: int = typer.Option(20, "--bam-min-baseq"),
    bam_min_mapq: int = typer.Option(20, "--bam-min-mapq"),
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
        policy = build_consensus_policy(
            bam_min_depth, bam_call_fract, bam_het_char, bam_lowcov, bam_min_baseq, bam_min_mapq
        )
        bam_resolver, _n_bam = load_bam_resolver(bam, bam_map)
        _run(sv_file, ref_genome, genome_build, flank, pop_vcf_dir, pop_data,
             pop_source, af_threshold, output, cols, bam_resolver, policy)
    except VflankError as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        raise typer.Exit(1) from exc


def _run(sv_file, ref_genome, genome_build, flank, pop_vcf_dir, pop_data,
         pop_source, af_threshold, output, cols: SvColumns, bam_resolver, policy):
    t0 = time.time()
    console.rule("[bold blue]vflank fusion run[/bold blue]")
    echo_parameters({
        "vflank version": __version__,
        "Breakpoints": sv_file, "Reference": ref_genome, "Genome build": genome_build,
        "Flank": f"{flank} bp/partner", "AF threshold": af_threshold,
        "Masking": (f"{pop_source} ({pop_data})" if (pop_vcf_dir or pop_source == "api")
                    else "none"),
        "BAM consensus": (
            f"on (min-depth={policy.min_depth}, het={policy.het_char}, "
            f"low-cov={policy.lowcov})" if bam_resolver is not None else "off"
        ),
        "Output": output,
    })
    if genome_build not in ("hg19", "hg38"):
        raise VflankError(f"--genome-build must be 'hg19' or 'hg38', got '{genome_build}'")
    validate_pop_options(pop_source, pop_data)
    if pop_vcf_dir is not None and not pop_vcf_dir.is_dir():
        raise VflankError(f"--pop-vcf-dir is not a directory: {pop_vcf_dir}")

    console.print(f"[bold]Loading breakpoints:[/bold] {sv_file}")
    df = bp_io.load_sv_table(sv_file, cols)
    console.print(f"  {len(df):,} fusion(s)")

    reference = ReferenceFasta(ref_genome)
    console.print(f"[bold]Reference:[/bold] {ref_genome}  [dim]({genome_build})[/dim]")
    build_warn = reference.check_build(genome_build)
    if build_warn:
        console.print(f"  [bold yellow]⚠ {build_warn}[/bold yellow]")

    # --- Masking source (optional) — masks the junction flanks ---
    bp_chroms = {
        b
        for col in (cols.chr1, cols.chr2)
        for b, _err in (normalise_chrom(v) for v in df[col].dropna().unique())
        if b
    }
    gnomad = make_pop_source(pop_source, pop_vcf_dir, genome_build, pop_data, bp_chroms)
    if pop_source == "api":
        dataset = dataset_for_build(genome_build)[1]
        console.print(f"[bold]Masking:[/bold] gnomAD API  [dim]({pop_data}, {dataset})[/dim]")
    elif gnomad is not None:
        console.print(f"[bold]Masking:[/bold] {pop_vcf_dir}  [dim](pop-data={pop_data})[/dim]")

    # --- BAM consensus (per-fusion patient sequence) ---
    bam_mode = bam_resolver is not None
    bam_cache: dict[str, object] = {}   # sample -> BamConsensusSource | None
    bam_warned: set[str] = set()
    n_consensus = 0
    if bam_mode:
        console.print(
            f"[bold]BAM consensus:[/bold] on  [dim](min-depth={policy.min_depth}, "
            f"het={policy.het_char}, low-cov={policy.lowcov})[/dim]"
        )
    console.print(f"[bold]Flank:[/bold] {flank} bp/partner (junction ≤ {2 * flank} bp)\n")

    def _bam_for(sample: str):
        if not bam_mode or not sample:
            return None
        if sample in bam_cache:
            return bam_cache[sample]
        path = bam_resolver(sample)
        src = None
        if path is not None:
            try:
                src = BamConsensusSource(path, policy)
            except ConsensusError as exc:
                if sample not in bam_warned:
                    log.warning("BAM unusable for %s (%s) — gnomAD masking only", sample, exc)
                    bam_warned.add(sample)
        elif sample not in bam_warned:
            log.warning("No BAM for sample %s — gnomAD masking only", sample)
            bam_warned.add(sample)
        bam_cache[sample] = src
        return src

    records: list[str] = []
    skipped = 0
    n_masked_total = 0
    skip_reasons: list[str] = []
    summary_rows: list[dict] = []

    for row_idx, row in df.iterrows():
        fusion, reason = bp_io.parse_fusion_row(row, cols)
        if reason is not None:
            skip_reasons.append(f"row {row_idx} — {reason}")
            skipped += 1
            continue
        bam_source = _bam_for(fusion.sample)
        try:
            jr = build_junction(
                reference, fusion, flank, gnomad=gnomad,
                af_threshold=af_threshold, bam_source=bam_source,
            )
        except Exception as exc:  # noqa: BLE001
            skip_reasons.append(f"row {row_idx} {fusion.name} — junction error: {exc}")
            skipped += 1
            continue
        if bam_source is not None:
            n_consensus += 1

        label = fasta_io.safe_header(fusion.name or "fusion")
        bp = f"{fusion.bp1.chrom}_{fusion.bp1.pos}_{fusion.bp1.strand}__" \
             f"{fusion.bp2.chrom}_{fusion.bp2.pos}_{fusion.bp2.strand}"
        prefix = f"{fasta_io.safe_header(fusion.sample)}__" if fusion.sample else ""
        header = f"{prefix}{label}__{bp}__j{jr.junction_index}"
        records.append(f">{header}\n{jr.sequence}\n")
        records.append(f">Masked__{header}\n{jr.masked_sequence}\n")
        n_masked_total += jr.n_masked

        truncated = len(jr.sequence) < 2 * flank
        summary_rows.append({
            "Name": fusion.name or ".", "BP1": f"{fusion.bp1.chrom}:{fusion.bp1.pos}",
            "BP2": f"{fusion.bp2.chrom}:{fusion.bp2.pos}",
            "Len": len(jr.sequence), "Junction": jr.junction_index,
            "N": jr.n_masked, "Trunc": truncated,
        })

    reference.close()
    if gnomad is not None:
        gnomad.close()
    for src in bam_cache.values():
        if src is not None:
            src.close()
    fasta_io.write_fasta(output, records)

    console.rule("[bold green]Results[/bold green]")
    if summary_rows:
        table = Table(show_header=True, header_style="bold cyan")
        for col in ("Name", "BP1", "BP2", "Len", "Junction", "N", "Trunc"):
            table.add_column(col)
        for r in summary_rows[:50]:
            n_str = f"[yellow]{r['N']}[/yellow]" if r["N"] else "[dim]0[/dim]"
            table.add_row(
                r["Name"], r["BP1"], r["BP2"], str(r["Len"]),
                str(r["Junction"]), n_str, "[yellow]yes[/yellow]" if r["Trunc"] else "no",
            )
        console.print(table)

    if skip_reasons:
        console.print(f"\n[bold yellow]Skipped {skipped}:[/bold yellow]")
        for reason in skip_reasons[:20]:
            console.print(f"  • {reason}")

    mask_line = f"[bold]Bases masked:[/bold] {n_masked_total:>6,}\n" if n_masked_total else ""
    consensus_line = f"[bold]BAM consensus:[/bold] {n_consensus:>6,}\n" if bam_mode else ""
    console.print(
        f"\n[bold]Fusions:[/bold]   {len(df):>6,}\n"
        + consensus_line +
        f"[bold]Records:[/bold]   {len(records):>6,} [dim](raw + masked per fusion)[/dim]\n"
        f"[bold]Skipped:[/bold]   {skipped:>6,}\n"
        + mask_line +
        f"[bold]Output:[/bold] [cyan]{output.resolve()}[/cyan]  "
        f"[dim]({time.time() - t0:.1f}s)[/dim]"
    )
