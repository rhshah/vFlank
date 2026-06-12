# Small Variants

`vflank small` extracts flanking sequence around SNPs/indels from a MAF and
writes a FASTA suitable for ddPCR assay design.

## Input

A tab-separated **MAF** (TCGA/MSK style). Required columns:
`Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`,
`Tumor_Seq_Allele2`. Optional metadata used in headers: `Hugo_Symbol`,
`HGVSp_Short`, `HGVSc`. Column names can be remapped with `--chrom-col`,
`--start-col`, … if your file differs.

Chromosome notation (`chr7` vs `7`) is auto-detected and normalised, including
numeric (`23`→X) and float (`17.0`→`17`, common when a column has blanks) forms.

```bash
vflank small inspect variants.maf      # preview columns + flag missing fields
```

## Run

```bash
vflank small run variants.maf \
    --ref-genome GRCh37.fasta \
    --genome-build hg19 \
    --pop-vcf-dir gnomad_v2.1.1/ \   # optional: SNP masking (or --pop-source api)
    --flank 200 \
    --output flanking_sequences.fasta \
    --report run_report.tsv          # optional machine-readable summary
```

## Output

Two records per **unique** variant — raw and `Masked__` — with the variant shown
literally as `[REF/ALT]` between the flanks. The header is keyed on the variant
identity, **not** the sample:

```
>{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
```

### Deduplication

The same variant seen across multiple samples collapses to **one** record
(flank + mask are sample-independent for reference/population masking). Use
`--no-dedup` to emit one record per row instead. The run summary reports how
many duplicates were collapsed.

## Safety nets

- **Genome-build guard** — if the FASTA's chr1 length disagrees with
  `--genome-build`, vflank warns (catches hg19-vs-hg38 mix-ups).
- **Flank truncation** — flanks that run off a contig end are emitted but
  reported, never silently shortened.
- **Skip summary** — invalid/incomplete rows (e.g. a missing `Chromosome`) are
  skipped and grouped by reason, with examples and a full list in `--report`.

## Sample filtering

```bash
--samples "P-001,P-002"        # comma-separated barcodes
--samples-file ids.txt          # one ID per line (# comments allowed)
```

See [SNP Masking](masking.md) for `--pop-source` / `--pop-data`.
