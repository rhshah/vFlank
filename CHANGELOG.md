# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Primer3 emit** (`--emit-primer3 FILE`) for both `small` and `fusion`: writes
  Boulder-IO records (`SEQUENCE_ID`/`TEMPLATE`/`TARGET`, and
  `SEQUENCE_EXCLUDED_REGION` from the masked positions) so a probe/primer design
  avoids common SNPs and patient-ambiguous sites. The template uses the
  masked/consensus call, filling `N` sites with the reference base and excluding
  them — a hard constraint, unlike a degenerate base. The `version` is recorded
  in the run report and parameter echo.

## [0.2.0] - 2026-06-14

### Added
- **Patient consensus from a BAM** (`--bam` / `--bam-map`) for both `small` and
  `fusion`: the masked record is the per-sample consensus (hom-ALT corrected,
  het/low-coverage handled), via `samtools consensus`. With a BAM, output is one
  record per (variant, sample) with the sample in the header; low coverage falls
  back to reference + gnomAD masking.
- **Coverage controls** for BAM consensus: `--bam-min-depth`, `--bam-call-fract`,
  `--bam-het-char {N,iupac}`, `--bam-lowcov`, `--bam-min-baseq`, `--bam-min-mapq`,
  and `--require-coverage` to flag (rather than silently fall back on) variants
  below the depth threshold.
- **Patient insertions are flagged, not dropped.** The reference-length consensus
  now scans read CIGARs and masks insertion anchor sites to `N`, surfaced as
  `NInserted` per variant and an "Insertion sites" line in the run summary —
  closing a silent-failure gap (deletions were already masked).
- **Run transparency:** every run echoes the full parameter set, and the TSV
  report carries per-variant detail (masked / corrected / inserted / coverage
  fraction / source / flagged).
- **Documentation:** glossary + site-wide abbreviation tooltips, mermaid
  architecture/flow diagrams (with pan/zoom), a left-sidebar layout, a
  "Highlighter" slate+amber palette, and a "Reading a record" explainer that
  shows masked sequences with inline highlighting, before/after diffs, and
  clickable annotations.

### Changed
- Adopted **git-flow** (`main` = releases, `develop` = integration) and reworked
  CI accordingly: a docs-build guard on every PR, versioned docs deployed with
  `mike` (`dev` from `develop`, `latest` from tags), and the GitHub Actions
  bumped off the deprecated Node-20 runtime (`checkout@v5`, `setup-python@v6`).

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

[Unreleased]: https://github.com/rhshah/vFlank/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/rhshah/vFlank/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/rhshah/vFlank/releases/tag/v0.1.0
