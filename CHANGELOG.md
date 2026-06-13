# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Patient consensus from a BAM** (`--bam` / `--bam-map`) for both `small` and
  `fusion`: the masked record is the per-sample consensus (hom-ALT corrected,
  het/low-coverage handled), via `samtools consensus`. With a BAM, output is one
  record per (variant, sample) with the sample in the header; low coverage falls
  back to reference + gnomAD masking.

## [0.1.0] - 2026-06-12

### Added
- Small-variant flank extraction from MAF (`vflank small run`), with `inspect`
  and `list-vcf` helpers.
- Common-SNP masking from two interchangeable backends (`--pop-source`):
  local gnomAD VCFs and the gnomAD GraphQL API (no download).
- `--pop-data {genome,exome,both}` for both backends.
- Per-variant deduplication keyed on `CHR_POS_REF_ALT` (`--dedup/--no-dedup`).
- Structural-variant junction extraction (`vflank fusion run`) from the simple
  iCallSV/iAnnotateSV breakpoint TSV (columns matched by name), with
  reverse-complement-aware junction construction and optional flank masking.
- Genome-build guard (hg19/hg38 vs FASTA), flank-truncation detection, and a
  categorised skip summary + optional TSV run report.
- Documentation site (MkDocs Material) and GitHub Actions CI.

[Unreleased]: https://github.com/rhshah/vFlank/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/rhshah/vFlank/releases/tag/v0.1.0
