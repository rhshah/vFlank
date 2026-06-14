# Patient Consensus from a BAM

gnomAD masking is a *population* prior. If you have the patient's aligned reads,
vflank can build the flank from the **patient's own sequence** — so a primer/probe
matches the real template, and `N` is used only where the patient is genuinely
ambiguous. This catches **private / rare** variants gnomAD never saw.

Available for both `vflank small` and `vflank fusion` via the same `--bam` /
`--bam-map` and `--bam-*` policy flags. (`--require-coverage`, which *flags*
low-coverage variants instead of falling back, is `vflank small` only.)

## How it works

For each flank position, vflank calls `samtools consensus` over the patient's
reads and decides:

| Patient genotype (depth ≥ `--bam-min-depth`) | Output base |
|---|---|
| homozygous-reference | reference base |
| **homozygous-alt** | the **patient's base** (consensus correction) |
| heterozygous | `N` (or an IUPAC code) |
| **low coverage** (`< --bam-min-depth`) | reference + gnomAD masking (see below) |

The patient consensus replaces the `Masked__` record; the raw record stays the
reference for comparison.

## Usage

=== "Single sample"

    ```bash
    vflank small run variants.maf -r GRCh37.fasta -g hg19 --bam sample.bam
    ```

=== "Cohort"

    A TSV mapping `Tumor_Sample_Barcode` to a BAM path:

    ```
    P-0001-T01-IM6	/data/P-0001-T01-IM6.bam
    P-0002-T01-IM6	/data/P-0002-T01-IM6.bam
    ```

    ```bash
    vflank small run variants.maf -r GRCh37.fasta -g hg19 --bam-map samples.tsv
    # fusion: identical flags; uses the breakpoint TSV's `sample` column
    vflank fusion run breakpoints.tsv -r GRCh37.fasta -g hg19 --bam-map samples.tsv
    ```

!!! note "Output granularity changes with a BAM"
    Because the consensus is **patient-specific**, the same variant in different
    samples is *not* collapsed: dedup keys on `(variant, sample)` and the
    **sample is added to the FASTA header**
    (`>{SAMPLE}__{GENE}__…__{CHROM}_{POS}_{REF}_{ALT}`). A sample without a BAM
    falls back to reference + gnomAD masking (with a warning).

## Low coverage

Below `--bam-min-depth` (default **20×**) we can't trust the reads, so vflank
falls back to **reference + gnomAD masking** by default — an uncovered variant
just behaves like a normal no-BAM run (it is *not* blanked to `N`). Combine with
`--pop-source` to N-mask common SNPs in those low-coverage stretches:

```bash
vflank small run variants.maf -r GRCh37.fasta -g hg19 \
    --bam-map samples.tsv --pop-source api      # consensus where covered, REF+gnomAD where not
```

`--bam-lowcov` overrides this: `n` (mask), `reference`, or `gnomad` (default).

## Options

| Option | Default | Meaning |
|---|---|---|
| `--bam` / `--bam-map` | — | single BAM / sample→BAM TSV |
| `--bam-min-depth` | 20 | minimum depth to trust a base |
| `--bam-call-fract` | 0.9 | fraction of reads to call a homozygous base |
| `--bam-het-char` | N | het output: `N` or `iupac` |
| `--bam-lowcov` | gnomad | low-coverage base: `n` / `reference` / `gnomad` |
| `--bam-min-baseq` / `--bam-min-mapq` | 20 / 20 | quality filters |

The BAM must be aligned to the same `--genome-build` as the reference and be
indexed (`.bai`). The engine is `samtools consensus` (bundled with pysam; no
external install). For large cohorts, parallelise per-sample or per-region with
Nextflow rather than in one process.
