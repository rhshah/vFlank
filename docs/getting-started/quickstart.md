# Quick Start

## Small variants (MAF → flanks)

```bash
vflank small run variants.maf \
    --ref-genome GRCh37.fasta \
    --genome-build hg19 \
    --flank 200 \
    --output flanking_sequences.fasta
```

Each unique variant yields two FASTA records — raw and masked:

```
>{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
{left_flank}[REF/ALT]{right_flank}
>Masked__{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
{left_flank_masked}[REF/ALT]{right_flank_masked}
```

Add SNP masking with **no download** via the gnomAD API:

```bash
vflank small run variants.maf -r GRCh37.fasta -g hg19 --pop-source api
```

Or skip the local FASTA too — fetch the reference from the UCSC API with
`--ref-source api` (then `--ref-genome` is not needed):

```bash
vflank small run variants.maf -g hg19 --ref-source api --pop-source api
```

!!! tip "Preview your MAF first"
    `vflank small inspect variants.maf` shows the columns and flags missing
    required fields before a full run.

## Fusions (breakpoint table → junctions)

Input is a tab-separated breakpoint table — columns matched **by name**:

```
chr1	pos1	str1	chr2	pos2	str2	name
22	29684066	0	11	32416521	1	EWSR1-WT1
```

(`str`: `0` = plus/reference, `1` = minus/complement, matching iCallSV.)

```bash
vflank fusion run breakpoints.tsv \
    --ref-genome GRCh37.fasta \
    --genome-build hg19 \
    --flank 200 \
    --output fusion_junctions.fasta
```

Each fusion yields the junction sequence (raw + masked), with the probe designed
to span the join.

## Next steps

- [Small Variants](../user-guide/small-variants.md)
- [Fusions / SVs](../user-guide/fusions.md)
- [SNP Masking](../user-guide/masking.md)
