# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.6.0] - 2026-06-17

### Added
- **VCF/BCF input for small variants** (`vflank small run variants.vcf.gz`). The
  input format is auto-detected by extension — no new command or flag. Read
  sites-only; anchor-base `REF`/`ALT` is normalised to MAF `[Start, End]`;
  multi-allelic records expand per `ALT`; symbolic/SV/BND alleles are skipped;
  gene/HGVS come best-effort from a VEP `CSQ` / SnpEff `ANN` field. `--samples`
  doesn't apply to a (sites-only) VCF and is ignored with a warning;
  `vflank small inspect` previews the normalised variants. New `io/vcf.py`
  (`is_vcf_path`, `vcf_to_maf_coords`, `load_vcf`); `run_small` accepts a VCF
  path too. SV-VCF (Delly `CT` / Manta–GRIDSS `BND`) remains a follow-up.

## [0.5.0] - 2026-06-15

### Added
- **Library API.** `vflank.pipeline.run_small` / `run_fusion(input, *, …) ->
  RunResult` run the whole pipeline (build sources → load → orchestrate → close)
  and return records, per-variant rows, categorised skips, counts, and API
  request tallies — no printing, no files written. Plus the streaming primitives
  `iter_small` / `iter_fusion` + `collect` for incremental use. This is the
  groundwork for hosting vflank as a service.
- **Buffer input.** `load_maf` / `load_sv_table` accept a path *or* an open
  text/binary buffer, so a service needn't round-trip through a temp file.

### Changed (internal)
- Extracted the per-variant orchestration out of the CLI into a presentation-free
  `vflank.pipeline` use-case layer (no Typer/Rich); the CLIs are now thin shells
  over it. Source construction moved to `vflank.sources`. CLI behaviour is
  unchanged. `pipeline.py` is now type-checked in the gate.

## [0.4.0] - 2026-06-14

### Added
- **No-download reference via the UCSC API** (`--ref-source api`) for both
  `small` and `fusion`: fetch each flank window from
  `api.genome.ucsc.edu/getData/sequence` instead of a local FASTA, so a run needs
  no multi-GB reference on disk (the prerequisite for hosting vflank as a
  service). `--ref-source` defaults to `file`; `--ref-genome` is now optional and
  required only for the file backend (clear error otherwise — no silent
  fallback). The UCSC backend uses 0-based half-open coordinates (identical to
  the local path, so flank math is unchanged), throttles to ~1 request/second,
  caches per window, retries transient failures, and surfaces a reference-API
  request count in the run summary and `--report` stats. See
  `docs/research/genome-api.md`.

### Changed
- `--help` for `small`/`fusion` now groups options into labelled panels
  (input, reference, masking, BAM, output) for readability.

### Docs
- Added a **Recommended Usage** page, a web-app/hosting design note
  (`docs/research/web-app-and-hosting.md`), and the reference-API note
  (`docs/research/genome-api.md`); recorded the pluggable reference source in the
  architecture and developer guides; added PyPI/GHCR/DeepWiki badges to the
  README and docs landing page.

## [0.3.0] - 2026-06-14

### Added
- **Primer3 emit** (`--emit-primer3 FILE`) for both `small` and `fusion`: writes
  Boulder-IO records (`SEQUENCE_ID`/`TEMPLATE`/`TARGET`, and
  `SEQUENCE_EXCLUDED_REGION` from the masked positions) so a probe/primer design
  avoids common SNPs and patient-ambiguous sites. The template uses the
  masked/consensus call, filling `N` sites with the reference base and excluding
  them — a hard constraint, unlike a degenerate base.
- **Version provenance**: the vflank version is recorded in the run report
  (`# vflank_version`) and the parameter echo, and a root `--version` flag is
  available alongside the `version` subcommand.

### Docs
- Refreshed the architecture roadmap and module maps to the shipped state, added
  a glossary + abbreviation tooltips, mermaid diagrams (pan/zoom), a slate+amber
  palette, a "Reading a record" explainer, and the Olivar/Primer3 emit-formats
  research note.

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

[Unreleased]: https://github.com/rhshah/vFlank/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/rhshah/vFlank/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/rhshah/vFlank/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/rhshah/vFlank/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/rhshah/vFlank/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/rhshah/vFlank/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/rhshah/vFlank/releases/tag/v0.1.0
