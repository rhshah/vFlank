# Installation

## Requirements

- **Python ≥ 3.10**
- Linux or macOS (pysam wheels are not published for native Windows — use WSL or
  conda there).
- An **indexed reference FASTA** (`.fai` alongside it; `samtools faidx ref.fasta`).
- Optionally, gnomAD VCFs **or** internet access to the gnomAD API (for masking).

## From GitHub

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

- a reference FASTA matching your `--genome-build` (`hg19`/GRCh37 or `hg38`/GRCh38);
- for masking, either a directory of gnomAD per-chromosome VCFs (`--pop-source vcf`)
  or nothing at all when using `--pop-source api`.

See [SNP Masking](../user-guide/masking.md) for details and download commands.
