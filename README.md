# vflank

[![CI](https://github.com/rhshah/vFlank/actions/workflows/ci.yml/badge.svg)](https://github.com/rhshah/vFlank/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs--material-blue)](https://rhshah.github.io/vFlank/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-green.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/rhshah/vFlank)

**Variant-aware flanking-sequence extraction and masking for ddPCR assay design.**

`vflank` is the *front-end* of a ddPCR assay-design pipeline. It takes genomic
variants — small variants (SNPs/indels) and structural variants (fusions) — and
emits the sequence an assay is designed around: the masked flanks of each variant
or the chimeric junction of a fusion. Primer/probe design itself is delegated
downstream to established tools.

📖 **Documentation: <https://rhshah.github.io/vFlank/>**

## Features

- **Small variants** (`vflank small`) — ±N bp flanks from a MAF, raw + masked
  FASTA, deduplicated per unique variant (`CHR_POS_REF_ALT`).
- **Fusions / SVs** (`vflank fusion`) — reverse-complement-aware junction
  sequences from an iCallSV / iAnnotateSV breakpoint table (columns by name).
- **SNP masking, two backends** — local gnomAD VCFs *or* the gnomAD GraphQL API
  (no download), each with `--pop-data {genome,exome,both}`.
- **Patient consensus from a BAM** (`--bam`/`--bam-map`) — build the flank/junction
  from the patient's own reads (hom-ALT corrected, het/low-cov handled) so primers
  match the real template; for both small variants and fusions.
- **No silent failures** — genome-build guard, flank-truncation detection, and a
  categorised skip summary + optional TSV report.

Planned: VCF input (small + BND SV) and downstream emit formats.
See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).

## Install

```bash
pip install vflank                                   # from PyPI (released versions)
pip install git+https://github.com/rhshah/vFlank.git # latest from GitHub
# development:
git clone https://github.com/rhshah/vFlank.git && cd vFlank
pip install -e ".[dev]"
```

Requires Python ≥ 3.10 (Linux/macOS) and `pysam`, `pandas`, `typer`, `rich`.

### Docker

Images are published to GHCR on each release:

```bash
docker run --rm -v "$PWD:/data" ghcr.io/rhshah/vflank \
    small run /data/variants.maf -r /data/GRCh37.fasta -g hg19 -o /data/out.fasta
```

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

### Masking sources

Common-SNP masking can come from local gnomAD VCFs or the gnomAD API:

- `--pop-source vcf` (default) — local per-chromosome gnomAD VCFs in
  `--pop-vcf-dir`. Reproducible, offline, unlimited scale.
- `--pop-source api` — the public [gnomAD GraphQL API](https://gnomad.broadinstitute.org/api),
  **no download**. Best for small cohorts (rate-limited to ~10 requests/min).

```bash
# No-download masking via the API (small cohorts):
vflank small run variants.maf -r GRCh37.fasta -g hg19 --pop-source api
```

Either source honours `--pop-data {genome,exome,both}` (default `genome`).
`both` masks a position if it is a common SNP in *either* the genome or exome
cohort. Flanks often fall in non-coding regions where only genomes have data,
so `genome` is the default.

Each variant yields two FASTA records (the `__{CHROM}_{POS}_{REF}_{ALT}` suffix
is what keys deduplication; the `{SAMPLE}__` prefix appears only with `--bam`):

```
>[{SAMPLE}__]{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
{left_flank}[REF/ALT]{right_flank}
>Masked__[{SAMPLE}__]{GENE}__{HGVSp}__{HGVSc}__{CHROM}_{POS}_{REF}_{ALT}
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
