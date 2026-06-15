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

from .. import __version__
from ..core.chrom import detect_series_chr_style, normalise_chrom
from ..core.consensus import BamConsensusSource, ConsensusFlankSource
from ..core.flanks import ReferenceFlankSource
from ..core.popfreq import build_chrom_vcf_map, kinds_for
from ..core.popfreq_api import dataset_for_build
from ..core.skips import categorize_skip
from ..errors import ConsensusError, VflankError
from ..io import emit_primer3 as primer3_io
from ..io import fasta as fasta_io
from ..io import maf as maf_io
from ..io import report as report_io
from ..io.maf import MAF_CHR, MAF_SAMPLE, REQUIRED_MAF_COLS, MafColumns
from ..logging import console, get_logger
from ._bam import build_consensus_policy, load_bam_resolver
from ._masking import make_pop_source, validate_pop_options
from ._reference import make_reference_source, validate_ref_source
from ._ui import (
    PANEL_ADVANCED,
    PANEL_BAM,
    PANEL_FILTER,
    PANEL_MAF_COLS,
    PANEL_MASKING,
    PANEL_REFERENCE,
    echo_parameters,
)

app = typer.Typer(no_args_is_help=True)
log = get_logger()


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
    ref_genome: Path | None = typer.Option(
        None, "--ref-genome", "-r",
        help="Indexed reference FASTA (.fai required). Required unless --ref-source api.",
        rich_help_panel=PANEL_REFERENCE,
    ),
    ref_source: str = typer.Option(
        "file", "--ref-source",
        help="Reference backend: file (local FASTA, default) or api (UCSC, no download).",
        rich_help_panel=PANEL_REFERENCE,
    ),
    pop_vcf_dir: Path | None = typer.Option(
        None, "--pop-vcf-dir", "-d",
        help="Directory of per-chromosome gnomAD VCF bgz files. Omit to skip masking.",
        rich_help_panel=PANEL_MASKING,
    ),
    genome_build: str = typer.Option(
        "hg19", "--genome-build", "-g",
        help="hg19 (GRCh37, gnomAD v2.1.1) or hg38 (GRCh38, gnomAD v4).",
    ),
    flank: int = typer.Option(
        200, "--flank", "-f", min=1, max=10_000, help="Bases on each side of the variant."
    ),
    af_threshold: float = typer.Option(
        0.001, "--af-threshold", min=0.0, max=1.0, help="Min population AF to mask a SNP.",
        rich_help_panel=PANEL_MASKING,
    ),
    pop_data: str = typer.Option(
        "genome", "--pop-data",
        help="gnomAD data to mask against: genome (default), exome, or both (union).",
        rich_help_panel=PANEL_MASKING,
    ),
    pop_source: str = typer.Option(
        "vcf", "--pop-source",
        help="Masking backend: vcf (local gnomAD VCFs) or api (gnomAD GraphQL, no download).",
        rich_help_panel=PANEL_MASKING,
    ),
    output: Path = typer.Option(
        Path("flanking_sequences.fasta"), "--output", "-o", help="Output FASTA file."
    ),
    report: Path | None = typer.Option(
        None, "--report",
        help="Write a per-variant TSV run report (stats + table) to this path.",
    ),
    emit_primer3: Path | None = typer.Option(
        None, "--emit-primer3",
        help="Also write a Primer3 Boulder-IO input file (one record per variant).",
    ),
    samples: str | None = typer.Option(
        None, "--samples", "-s",
        help="Comma-separated Tumor_Sample_Barcode IDs to include.",
        rich_help_panel=PANEL_FILTER,
    ),
    samples_file: Path | None = typer.Option(
        None, "--samples-file",
        help="File of sample IDs, one per line (# comments allowed).",
        rich_help_panel=PANEL_FILTER,
    ),
    chrom_col: str = typer.Option(MafColumns.chrom, "--chrom-col", rich_help_panel=PANEL_MAF_COLS),
    start_col: str = typer.Option(MafColumns.start, "--start-col", rich_help_panel=PANEL_MAF_COLS),
    end_col: str = typer.Option(MafColumns.end, "--end-col", rich_help_panel=PANEL_MAF_COLS),
    ref_col: str = typer.Option(MafColumns.ref, "--ref-col", rich_help_panel=PANEL_MAF_COLS),
    alt_col: str = typer.Option(MafColumns.alt, "--alt-col", rich_help_panel=PANEL_MAF_COLS),
    gene_col: str = typer.Option(MafColumns.gene, "--gene-col", rich_help_panel=PANEL_MAF_COLS),
    prot_col: str = typer.Option(MafColumns.protein, "--prot-col", rich_help_panel=PANEL_MAF_COLS),
    cdna_col: str = typer.Option(MafColumns.cdna, "--cdna-col", rich_help_panel=PANEL_MAF_COLS),
    sample_col: str = typer.Option(
        MafColumns.sample, "--sample-col", rich_help_panel=PANEL_MAF_COLS
    ),
    uppercase: bool = typer.Option(
        True, "--uppercase/--no-uppercase", help="Uppercase flanking sequences.",
        rich_help_panel=PANEL_ADVANCED,
    ),
    dedup: bool = typer.Option(
        True, "--dedup/--no-dedup",
        help="Emit one record per unique variant (CHR_POS_REF_ALT), collapsing samples.",
        rich_help_panel=PANEL_ADVANCED,
    ),
    bam: Path | None = typer.Option(
        None, "--bam", help="Single-sample BAM for patient consensus (modes C/D).",
        rich_help_panel=PANEL_BAM,
    ),
    bam_map: Path | None = typer.Option(
        None, "--bam-map",
        help="TSV (Tumor_Sample_Barcode<TAB>bam_path) for cohort consensus.",
        rich_help_panel=PANEL_BAM,
    ),
    bam_min_depth: int = typer.Option(
        20, "--bam-min-depth", help="Min depth to trust a base.", rich_help_panel=PANEL_BAM
    ),
    bam_call_fract: float = typer.Option(
        0.9, "--bam-call-fract", help="Fraction to call a base.", rich_help_panel=PANEL_BAM
    ),
    bam_het_char: str = typer.Option(
        "N", "--bam-het-char", help="Het output: N or iupac.", rich_help_panel=PANEL_BAM
    ),
    bam_lowcov: str = typer.Option(
        "gnomad", "--bam-lowcov",
        help="Low-coverage base: n (mask) | reference | gnomad (default: REF + gnomAD).",
        rich_help_panel=PANEL_BAM,
    ),
    bam_min_baseq: int = typer.Option(20, "--bam-min-baseq", rich_help_panel=PANEL_BAM),
    bam_min_mapq: int = typer.Option(20, "--bam-min-mapq", rich_help_panel=PANEL_BAM),
    require_coverage: float = typer.Option(
        0.0, "--require-coverage", min=0.0, max=1.0,
        help="Flag BAM-consensus variants whose flanks are < this fraction covered (0=off).",
        rich_help_panel=PANEL_BAM,
    ),
):
    """Extract flanking sequences for every variant in a MAF and write a FASTA.

    Each variant yields a raw record and a Masked record. Common population SNPs
    (gnomAD, AF >= --af-threshold) are masked; with --bam/--bam-map the Masked
    record is the per-sample patient consensus instead (one record per
    (variant, sample); sample is added to the header). Chromosome notation is
    auto-detected from the FASTA and VCFs.
    """
    try:
        policy = build_consensus_policy(
            bam_min_depth, bam_call_fract, bam_het_char, bam_lowcov, bam_min_baseq, bam_min_mapq
        )
        bam_resolver, n_bam = load_bam_resolver(bam, bam_map)
        _run(
            maf_file, ref_genome, ref_source, pop_vcf_dir, genome_build, flank, af_threshold,
            pop_data, pop_source, output, report, emit_primer3, samples, samples_file,
            MafColumns(chrom_col, start_col, end_col, ref_col, alt_col,
                       gene_col, prot_col, cdna_col, sample_col),
            uppercase, dedup, bam_resolver, n_bam, policy, require_coverage,
        )
    except VflankError as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        raise typer.Exit(1) from exc


def _run(maf_file, ref_genome, ref_source, pop_vcf_dir, genome_build, flank, af_threshold,
         pop_data, pop_source, output, report, emit_primer3, samples, samples_file,
         cols: MafColumns, uppercase: bool, dedup: bool,
         bam_resolver, n_bam, policy, require_coverage):
    t0 = time.time()
    console.rule("[bold blue]vflank small run[/bold blue]")

    bam_mode = bam_resolver is not None
    echo_parameters({
        "vflank version": __version__,
        "MAF": maf_file,
        "Reference": (ref_genome if ref_source == "file" else "UCSC API"),
        "Genome build": genome_build,
        "Flank": f"±{flank} bp", "AF threshold": af_threshold,
        "Masking": (f"{pop_source} ({pop_data})" if (pop_vcf_dir or pop_source == "api")
                    else "none"),
        "Dedup": dedup, "Uppercase": uppercase,
        "Sample filter": samples or samples_file or "none",
        "BAM consensus": (
            f"on (min-depth={policy.min_depth}, call-fract={policy.call_fract}, "
            f"het={policy.het_char}, low-cov={policy.lowcov}, "
            f"baseq={policy.min_baseq}, mapq={policy.min_mapq})" if bam_mode else "off"
        ),
        "Require coverage": require_coverage or "off",
        "Output": output, "Report": report or "none",
        "Emit Primer3": emit_primer3 or "off",
    })

    if genome_build not in ("hg19", "hg38"):
        raise VflankError(f"--genome-build must be 'hg19' or 'hg38', got '{genome_build}'")
    validate_ref_source(ref_source)
    validate_pop_options(pop_source, pop_data)
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

    # --- Reference source (local FASTA or UCSC API) + build guard ---
    reference = make_reference_source(ref_source, ref_genome, genome_build)
    ref_label = "UCSC API" if ref_source == "api" else ref_genome
    console.print(f"[bold]Reference:[/bold] {ref_label}  [dim]({genome_build})[/dim]")
    build_warn = reference.check_build(genome_build)
    if build_warn:
        console.print(f"  [bold yellow]⚠ {build_warn}[/bold yellow]")

    maf_style = detect_series_chr_style(df[MAF_CHR].head(20)) if MAF_CHR in df.columns else None
    console.print(
        f"  [dim]FASTA contigs: {'chr-prefixed' if reference.has_chr else 'bare'} · "
        f"MAF: {'chr-prefixed' if maf_style else 'bare' if maf_style is False else 'unknown'}[/dim]"
    )

    # --- Masking source (optional) ---
    maf_chroms = {
        b for b, _err in (normalise_chrom(c) for c in df[MAF_CHR].dropna().unique()) if b
    }
    gnomad = make_pop_source(pop_source, pop_vcf_dir, genome_build, pop_data, maf_chroms)
    if pop_source == "api":
        if pop_vcf_dir is not None:
            console.print("  [yellow]⚠ --pop-vcf-dir ignored with --pop-source api[/yellow]")
        if len(df) > 50:
            console.print(
                f"  [yellow]⚠ {len(df):,} variants over the rate-limited API "
                f"(~10 req/min) — consider --pop-source vcf for bulk.[/yellow]"
            )
        dataset = dataset_for_build(genome_build)[1]
        console.print(
            f"[bold]Masking:[/bold] gnomAD API  "
            f"[dim](pop-data={pop_data}, dataset={dataset}, AF ≥ {af_threshold})[/dim]"
        )
    elif gnomad is not None:
        console.print(
            f"[bold]Population VCF dir:[/bold] {pop_vcf_dir}  "
            f"[dim](pop-data={pop_data}, AF ≥ {af_threshold})[/dim]"
        )
    else:
        console.print(
            "[yellow]⚠ No masking source (--pop-source vcf without --pop-vcf-dir).[/yellow]"
        )
    console.print(f"[bold]Flank:[/bold] ±{flank} bp")

    ref_flank_source = ReferenceFlankSource(
        reference, gnomad, flank=flank, af_threshold=af_threshold
    )

    # --- BAM consensus mode (per-sample patient sequence) ---
    # Default low-coverage behaviour is REF + gnomAD masking: where the BAM is
    # shallow (< min_depth) we fall back to the reference base, with gnomAD
    # N-masking common SNPs if a population source is given. So an uncovered
    # variant just behaves like a normal no-BAM run (not all-N).
    consensus_cache: dict[str, object] = {}   # sample -> ConsensusFlankSource (or ref fallback)
    bam_warned: set[str] = set()
    n_consensus = 0
    n_flagged = 0
    n_inserted_total = 0
    if bam_mode:
        scope = "single BAM (all samples)" if n_bam == -1 else f"{n_bam} sample(s) mapped"
        console.print(
            f"[bold]BAM consensus:[/bold] {scope}  [dim](min-depth={policy.min_depth}, "
            f"het={policy.het_char}, low-cov={policy.lowcov})[/dim]"
        )
    console.print()

    def _source_for(variant):
        """Pick the flank source for a variant (consensus per-sample, or reference)."""
        if not bam_mode:
            return ref_flank_source, False
        sample = variant.sample
        if sample in consensus_cache:
            cached = consensus_cache[sample]
            return cached, not isinstance(cached, ReferenceFlankSource)
        bam_path = bam_resolver(sample)
        if bam_path is None:
            if sample not in bam_warned:
                log.warning("No BAM for sample %s — using reference + gnomAD masking", sample)
                bam_warned.add(sample)
            consensus_cache[sample] = ref_flank_source
            return ref_flank_source, False
        try:
            src = ConsensusFlankSource(
                reference, BamConsensusSource(bam_path, policy), gnomad,
                flank=flank, af_threshold=af_threshold,
            )
        except ConsensusError as exc:
            if sample not in bam_warned:
                log.warning("BAM unusable for %s (%s) — reference + gnomAD", sample, exc)
                bam_warned.add(sample)
            consensus_cache[sample] = ref_flank_source
            return ref_flank_source, False
        consensus_cache[sample] = src
        return src, True

    # --- Process variants ---
    records: list[str] = []
    primer3_records: list[primer3_io.Primer3Record] = []
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

            # Dedup: by CHR_POS_REF_ALT (sample-independent reference/gnomAD
            # masking). With a BAM the consensus is patient-specific, so the key
            # also includes the sample -> one record per (variant, sample).
            if dedup:
                key = (variant.chrom, variant.start, variant.end,
                       variant.ref.upper(), variant.alt.upper())
                if bam_mode:
                    key = (*key, variant.sample)
                if key in seen_variants:
                    n_duplicate += 1
                    continue
                seen_variants.add(key)

            source, used_consensus = _source_for(variant)
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

            if used_consensus:
                n_consensus += 1
            n_masked_total += fr.n_masked
            sample_tag = variant.sample if bam_mode else None
            records.extend(fasta_io.format_records(variant, fr, ref, alt, sample=sample_tag))
            if emit_primer3 is not None:
                primer3_records.append(primer3_io.small_variant_record(
                    fasta_io.record_id(variant, ref, alt, sample_tag),
                    fr.left, fr.right, fr.masked_left, fr.masked_right, ref,
                ))

            row_detail = {
                "Sample": variant.sample[:32], "Gene": variant.gene,
                "Chrom": variant.raw_chrom, "Start": variant.start, "End": variant.end,
                "Ref": ref[:8], "Alt": alt[:8],
                "LeftLen": len(fr.left), "RightLen": len(fr.right),
                "NMasked": fr.n_masked, "NCorrected": fr.n_corrected, "Truncated": truncated,
            }
            if bam_mode:
                source_label = "consensus" if used_consensus else "reference"
                covered_frac = (fr.covered / fr.total) if (used_consensus and fr.total) else None
                flagged = (
                    require_coverage > 0 and covered_frac is not None
                    and covered_frac < require_coverage
                )
                if flagged:
                    n_flagged += 1
                row_detail["Source"] = source_label
                row_detail["CoveredFrac"] = (
                    round(covered_frac, 3) if covered_frac is not None else ""
                )
                row_detail["NInserted"] = fr.inserted or 0
                row_detail["Flagged"] = flagged
                n_inserted_total += fr.inserted or 0
            summary_rows.append(row_detail)

    ref_api_requests = getattr(reference, "request_count", None) if ref_source == "api" else None
    reference.close()
    api_requests = getattr(gnomad, "request_count", None) if gnomad is not None else None
    if gnomad is not None:
        gnomad.close()
    for src in consensus_cache.values():
        if hasattr(src, "bam"):
            src.bam.close()

    fasta_io.write_fasta(output, records)
    if emit_primer3 is not None:
        primer3_io.write_primer3(emit_primer3, primer3_records)

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
        f"[bold]gnomAD API req:[/bold]   {api_requests:>6,}\n" if api_requests is not None else ""
    )
    ref_api_line = (
        f"[bold]Reference API req:[/bold] {ref_api_requests:>5,}\n"
        if ref_api_requests is not None else ""
    )
    dup_line = (
        f"[bold]Dup. collapsed:[/bold]   {n_duplicate:>6,} [dim](other samples)[/dim]\n"
        if n_duplicate else ""
    )
    consensus_line = (
        f"[bold]BAM consensus:[/bold]    {n_consensus:>6,} [dim](patient-specific)[/dim]\n"
        if bam_mode else ""
    )
    inserted_line = (
        f"[bold]Insertion sites:[/bold]  {n_inserted_total:>6,} [dim](masked N)[/dim]\n"
        if (bam_mode and n_inserted_total) else ""
    )
    flagged_line = (
        f"[bold yellow]Low-coverage flagged:[/bold yellow] {n_flagged:>6,} "
        f"[dim](< {require_coverage:.0%} covered)[/dim]\n"
        if (bam_mode and require_coverage > 0) else ""
    )
    console.print(
        f"\n[bold]Total in MAF:[/bold]     {n_total:>6,}\n"
        f"[bold]Processed:[/bold]        {len(summary_rows):>6,}\n"
        f"[bold]Skipped:[/bold]          {skipped:>6,}\n"
        + dup_line + truncated_line + consensus_line + inserted_line + flagged_line +
        f"[bold]Bases masked:[/bold]     {n_masked_total:>6,}\n"
        + api_line + ref_api_line +
        f"[bold]FASTA records:[/bold]    {len(records):>6,} [dim](2 per variant)[/dim]\n"
        f"[bold]Output:[/bold] [cyan]{output.resolve()}[/cyan]  [dim]({elapsed:.1f}s)[/dim]"
    )

    # --- Optional machine-readable report ---
    if report is not None:
        stats = {
            # provenance
            "vflank_version": __version__,
            # run parameters (what was set)
            "maf": maf_file,
            "reference": (str(ref_genome) if ref_source == "file" else "UCSC API"),
            "ref_source": ref_source, "genome_build": genome_build,
            "flank": flank, "af_threshold": af_threshold,
            "pop_source": pop_source, "pop_data": pop_data, "dedup": dedup,
            # outcomes (what happened)
            "total_in_maf": n_total,
            "processed": len(summary_rows),
            "skipped": skipped,
            "duplicates_collapsed": n_duplicate,
            "truncated_flanks": n_truncated,
            "bases_masked": n_masked_total,
            "fasta_records": len(records),
        }
        if emit_primer3 is not None:
            stats["primer3_records"] = len(primer3_records)
        if api_requests is not None:
            stats["api_requests"] = api_requests
        if ref_api_requests is not None:
            stats["reference_api_requests"] = ref_api_requests
        if bam_mode:
            stats["bam_consensus_records"] = n_consensus
            stats["bam_min_depth"] = policy.min_depth
            stats["bam_het_char"] = policy.het_char
            stats["bam_lowcov"] = policy.lowcov
            stats["require_coverage"] = require_coverage
            stats["low_coverage_flagged"] = n_flagged
            stats["insertion_sites_masked"] = n_inserted_total
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
