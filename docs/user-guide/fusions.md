# Fusions / Structural Variants

`vflank fusion` builds the chimeric **junction sequence** for a fusion so a ddPCR
probe can span it.

## Input

A tab-separated **breakpoint table**, columns matched **by name** (any order,
extra columns ignored):

| Column | Meaning |
|--------|---------|
| `chr1`, `pos1`, `str1` | breakpoint 1: chromosome, 1-based position, strand |
| `chr2`, `pos2`, `str2` | breakpoint 2 |
| `name`, `sample` | optional labels |

`str` follows the iCallSV convention: **`0` = plus/reference, `1` = minus/complement**.
This is the format produced by [iCallSV](https://github.com/rhshah/iCallSV) /
[iAnnotateSV](https://github.com/rhshah/iAnnotateSV). Rename columns with
`--chr1-col`, `--pos1-col`, … if needed.

```
name	chr1	pos1	str1	chr2	pos2	str2
EWSR1-WT1	22	29684066	0	11	32416521	1
```

## Run

```bash
vflank fusion run breakpoints.tsv \
    --ref-genome GRCh37.fasta \
    --genome-build hg19 \
    --flank 200 \
    --pop-source api \           # optional: mask common SNPs in the junction flanks
    --output fusion_junctions.fasta
```

## How the junction is built

The fused product reads 5'→3' as `partner1 + partner2` (no separator). With
`L = --flank` bases per side:

| | partner 1 (ends at junction) | partner 2 (starts at junction) |
|---|---|---|
| **str = 0** | `ref[pos−L+1 .. pos]` (+) | `revcomp(ref[pos−L+1 .. pos])` |
| **str = 1** | `revcomp(ref[pos .. pos+L−1])` | `ref[pos .. pos+L−1]` (+) |

The probe is designed to span `junction_index` (the FASTA header records it as
`__j{index}`). Junctions that run off a contig end are emitted but flagged.

## Output

One raw + one `Masked__` record per fusion (the `{SAMPLE}__` prefix appears only
when a per-sample BAM consensus is used, mirroring `vflank small`):

```
>[{SAMPLE}__]{NAME}__{chr1}_{pos1}_{str1}__{chr2}_{pos2}_{str2}__j{index}
{junction}
```

## Masking

The same masking applies to the junction flanks (a real upgrade over per-fusion
scripts that ignore polymorphisms). Masking is computed in genomic space before
any reverse-complement — see [SNP Masking](masking.md).
