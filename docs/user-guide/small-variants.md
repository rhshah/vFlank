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

``` { .bash .annotate }
vflank small run variants.maf \
    --ref-genome GRCh37.fasta \
    --genome-build hg19 \            # (1)!
    --pop-vcf-dir gnomad_v2.1.1/ \   # (2)!
    --flank 200 \                    # (3)!
    --output flanking_sequences.fasta \
    --report run_report.tsv          # (4)!
```

1.  Guards against an hg19-vs-hg38 mix-up — vflank checks the FASTA's chr1
    length against this build and warns on mismatch.
2.  Optional SNP masking from local gnomAD VCFs. Swap for `--pop-source api` to
    use the gnomAD GraphQL API instead (no download).
3.  Bases taken from each side of the variant (so each record is up to
    `2 × flank` bp, shorter at a contig end).
4.  Optional machine-readable run summary: per-variant masked/corrected counts,
    skips grouped by reason, and the full parameter set.

## Output

Two records per **unique** variant — raw and `Masked__` — with the variant shown
literally as `[REF/ALT]` between the flanks. The header is keyed on the variant
identity, **not** the sample:

```
>{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
```

### Reading a record

A masked flank is just the reference sequence with the variant shown literally
and common SNPs swapped for `N`. Highlighting only the parts vflank touches (the
rest is untouched, designable reference):

AGCGATCGATCGTACGT==[T/C]==ACGTGCA==N==TCGATCGTAGC

…where ==[T/C]== is the variant of interest and ==N== marks a masked common SNP
(gnomAD AF ≥ threshold) that a primer/probe must avoid.

**Before → after masking** — what the `Masked__` record changes versus the raw
record (one common SNP, here, replaced):

```diff
- AGCGATCGATCGTACGTGCAATCGATCGTAGC
+ AGCGATCGATCGTACGTGCAATCGNTCGTAGC
```

**Anatomy** of a record, left to right:

`left flank (5′)`
:   `--flank` bases of reference ending just *before* the variant.

`[REF/ALT]`
:   the variant itself, written literally; excluded from both flanks.

`right flank (3′)`
:   `--flank` bases starting just *after* the variant — masked `N` wherever a
    common SNP would sit under a primer/probe.

### Deduplication

The same variant seen across multiple samples collapses to **one** record
(flank + mask are sample-independent for reference/population masking). Use
`--no-dedup` to emit one record per row instead. The run summary reports how
many duplicates were collapsed.

## Emit for Primer3

Add `--emit-primer3 primers.txt` to also write a [Primer3](https://primer3.org)
Boulder-IO input file — one record per variant, ready to hand to a designer:

```bash
vflank small run variants.maf -r GRCh37.fasta -g hg19 \
    -o flanks.fasta --emit-primer3 primers.txt
```

Each record carries:

- `SEQUENCE_TEMPLATE` — the best-known sequence (the masked/consensus call,
  falling back to the reference base where the call is `N`).
- `SEQUENCE_TARGET` — the variant span, so the assay covers it.
- `SEQUENCE_EXCLUDED_REGION` — the masked positions (common SNPs, patient
  het/low-cov/insertion sites). This is a **hard** "no oligo here" constraint —
  stronger than a degenerate `N`, which Primer3 may still design over.

`SEQUENCE_ID` matches the FASTA header key, so the two outputs cross-reference.
The same flag works on `vflank fusion run` (one record per junction, targeted to
span the breakpoint). Olivar emit is [planned](../research/emit-formats.md).

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
