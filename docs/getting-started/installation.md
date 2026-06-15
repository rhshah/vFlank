# Installation

## Requirements

- **Python ≥ 3.10**
- Linux or macOS (pysam wheels are not published for native Windows — use WSL or
  conda there).
- An **indexed reference FASTA** (`.fai` alongside it; `samtools faidx ref.fasta`)
  — *or* internet access to use `--ref-source api` (UCSC) with no local FASTA.
- Optionally, gnomAD VCFs **or** internet access to the gnomAD API (for masking).

## From PyPI

```bash
pip install vflank
```

## From GitHub (latest, unreleased)

```bash
pip install git+https://github.com/rhshah/vFlank.git
```

## From source (development)

```bash
git clone https://github.com/rhshah/vFlank.git
cd vFlank
pip install -e ".[dev]"
```

Runtime dependencies (`typer`, `rich`, `pysam`, `pandas`) install automatically.

## Verify

```bash
vflank version
vflank --help
```

You should see the `small` and `fusion` command groups.

## Reference and masking data

vflank never bundles large genomic data. You provide:

- a reference matching your `--genome-build` (`hg19`/GRCh37 or `hg38`/GRCh38):
  either a local indexed FASTA (`--ref-source file`, default) or the UCSC API
  (`--ref-source api`, no download);
- for masking, either a directory of gnomAD per-chromosome VCFs (`--pop-source vcf`)
  or nothing at all when using `--pop-source api`.

See [SNP Masking](../user-guide/masking.md) for details and download commands.
