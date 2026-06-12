# vflank

**Variant-aware flanking-sequence extraction and masking for ddPCR assay design.**

`vflank` is the *front-end* of a ddPCR assay-design pipeline. It takes genomic
variants, extracts the reference (and, soon, patient-specific) sequence flanking
each one, masks positions that would sabotage a primer/probe (common germline
SNPs, and — soon — patient-specific heterozygous sites observed in a BAM), and
emits clean target sequences. Primer/probe design itself is delegated downstream
to established tools (Olivar for small-variant amplicons, Primer3 for
fusion-junction probes).

## Status

Early development. Implemented today:

- `vflank small run` — extract ± *N* bp flanks for every variant in a MAF, mask
  common gnomAD SNPs (AF ≥ threshold), write a raw + masked FASTA per variant.
- `vflank small inspect` — preview MAF columns.
- `vflank small list-vcf` — verify gnomAD per-chromosome coverage.

Planned: BAM consensus flanks (modes C/D), fusion-junction rewrite (Python 3 +
Pydantic config), `--emit-olivar`/`--emit-primer3` outputs, and a Nextflow
pipeline wrapping the CLIs in containers. See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).

## Install (development)

```bash
pip install -e ".[dev]"     # needs the hatchling build backend
# or run tests without installing:
python -m pytest
```

Requires Python ≥ 3.10 and `pysam`, `pandas`, `typer`, `rich`.

## Quick start

```bash
vflank small run variants.maf \
    --ref-genome /path/to/GRCh37.fasta \
    --pop-vcf-dir /path/to/gnomad_v2.1.1/ \
    --genome-build hg19 \
    --flank 200 \
    --output flanking_sequences.fasta
```

`--genome-build` defaults to **hg19** (GRCh37 / gnomAD v2.1.1); pass `-g hg38`
for GRCh38 / gnomAD v4. gnomAD v4 has no GRCh37 build.

Each variant yields two FASTA records:

```
>{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}
{left_flank}[REF/ALT]{right_flank}
>Masked__{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}
{left_flank_masked}[REF/ALT]{right_flank_masked}
```

Chromosome notation (`chr1` vs `1`) is auto-detected from the FASTA and VCFs.
The genome build is sanity-checked against the FASTA's chr1 length.

## Project layout

```
src/vflank/
├── core/   chrom · variant · flanks · popfreq   (pure, testable domain logic)
├── io/     maf · reference · fasta              (file access)
└── cli/    app · small                          (Typer commands)
```

## Documentation

- [docs/DEVELOPER.md](docs/DEVELOPER.md) — setup, running, testing, using vflank
  as a library, and extending it (new flank sources, CLI commands).
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — design, scope boundary, and the
  milestone roadmap.
- `CLAUDE.md` — repository conventions and the quality gate.
