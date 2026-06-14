# CLAUDE.md — working guide for the vflank repository

## What this project is

`vflank` is the **variant-aware, masked-flank front-end** of a ddPCR
assay-design pipeline. It extracts the sequence flanking genomic variants and
masks positions that would compromise a primer/probe, then emits clean target
sequences.

It is **not** a primer designer. Design is delegated downstream — Olivar
(small-variant amplicons) and Primer3 (fusion-junction probes) — invoked
out-of-process. Do not add a primer/probe design algorithm to this package; add
*emit formats* that feed those tools. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)
for the full plan, scope boundary, and milestone roadmap.

## Repository layout

```
src/vflank/
├── core/   chrom · variant · flanks · popfreq · popfreq_api · fusion · skips
├── io/     maf · reference · fasta · breakpoints · report   ← file access
├── cli/    app (root) · small (run/inspect/list-vcf) · fusion (run)
├── logging.py · errors.py
tests/      unit/ · integration/
docs/ARCHITECTURE.md                              ← design & roadmap
docs/research/                                     ← gnomAD-API & SV/VCF design notes
```

The original reference scripts (`get_flanking_sequence.py`,
`design_fusion_primers.py`, `config_ES_CTDNA_03.cfg`) have been fully ported /
re-implemented and removed. The conventions extracted from them live in
[docs/research/](docs/research/) (corrected fusion-junction model, iCallSV
strand mapping, etc.).

## Quality gate — run before declaring any change done

```bash
python -m ruff check src tests
python -m mypy src/vflank/core src/vflank/io
python -m pytest
```

All three must pass. Tests run without installing the package (`pyproject.toml`
sets `pythonpath = ["src"]`). The dev environment is mambaforge Python 3.10 with
`typer`, `rich`, `pysam`, `pandas`, `pytest`, `ruff`, `mypy` available.

## Working discipline (the bar for this repo)

- **Review before and after every change.** Before: read the surrounding code
  and run the gate. After: re-run the gate; re-read the diff for duplication,
  dead code, and unused symbols.
- **No silent failures.** Every error path must surface — raise a typed
  `VflankError`, or log at an appropriate level, or record and report in the run
  summary. Do not swallow exceptions or return empty results without a log. The
  existing patterns: flank truncation at contig ends is detected and reported;
  contig-absent in a VCF is logged at DEBUG; build mismatch warns.
- **No dead or duplicated code.** If you remove the last caller of something,
  remove the thing. Don't add an exception class / helper "for later."
- **Keep the hot kernels pure.** `chrom.normalise_chrom`, `popfreq.parse_common_snp_positions`,
  and `flanks.mask_sequence` are pure functions over plain values so they are
  unit-testable without pysam and can later be ported to Rust. Preserve that.
- **Update comments, logging, and tests with the code** — not as an afterthought.

## Coordinate conventions (the #1 source of bugs — read before editing flanks)

- **MAF coordinates are 1-based, fully-closed `[start, end]`.**
- **pysam (`FastaFile.fetch`, tabix) is 0-based, half-open `[start, end)`.**
- Flank math (see `core/flanks.ReferenceFlankSource.fetch`):
  - left flank  = `fa.fetch(chrom, max(0, start-flank-1), start-1)` — bases *before* the variant
  - right flank = `fa.fetch(chrom, end, end+flank)` — bases *after* the variant
  - the variant interval `[start, end]` itself is excluded from both flanks
- Masking maps a 1-based VCF position back to a 0-based index *within a flank*:
  `idx = pos - region_start_0based - 1` (see `flanks.mask_sequence`).
- pysam silently returns a **short string** when a window runs off a contig end —
  always treat a shorter-than-requested flank as a real condition to report.

## Domain knowledge

For the deeper biology (why masking matters for ddPCR, gnomAD AF semantics across
builds, masking modes A–D, indel/het caveats), invoke the `ddpcr-conventions`
skill — it carries the reference detail that doesn't belong in always-on context.

## Conventions

- **Chromosome handling:** the canonical internal form is the *bare* chromosome
  (`"7"`, `"X"`, `"MT"`). Normalise MAF input with `core/chrom.normalise_chrom`;
  convert to a file's notation only at fetch time via `ReferenceFasta.contig` /
  the gnomAD store. Notation (`chr1` vs `1`) is auto-detected per file.
- **Output:** two FASTA records per variant — raw and `Masked__…` — with the
  variant shown literally as `[REF/ALT]` between the flanks. Only single-base
  population SNPs are masked.
- **CLI:** Typer; status/Rich output goes to **stderr** (the `console` in
  `logging.py`); stdout/files are for data. Library/diagnostic messages go
  through the `vflank` logger; the CLI presents the formatted summary.
- **Errors:** raise `VflankError` subclasses (`errors.py`) for user-facing
  failures; the CLI catches them and prints a clean message.

## Git

The repo follows **git-flow** (`main` = production, `develop` = integration):

- **`main`** holds released, tagged versions only (`vX.Y.Z`); never commit to it
  directly. **`develop`** is the integration branch — feature work branches off
  it and merges back into it.
- Day-to-day work uses git-flow branches: `git flow feature start <name>`
  (branches off `develop`), `git flow release start X.Y.Z`, `git flow hotfix
  start X.Y.Z` (off `main`). Prefixes: `feature/ release/ hotfix/`; version-tag
  prefix `v`. (Installed git-flow is nvie 0.4.1 — no `bugfix` subcommand; use a
  `feature` branch for fixes.)
- A release/hotfix `finish` merges into **both** `main` and `develop` and tags
  `main`. Small doc-only tweaks may be committed straight to `develop`.
- End commit messages with: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
- Commit/push only when asked.
