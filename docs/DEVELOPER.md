# vflank Developer Guide

How to set up, run, test, extend, and reuse `vflank`. For the *why* (design and
roadmap) see [ARCHITECTURE.md](ARCHITECTURE.md); for the biology see the
`ddpcr-conventions` skill; for the user-facing quickstart see the top-level
[README](https://github.com/rhshah/vFlank/blob/main/README.md).

---

## 1. Set up a development environment

```bash
# From the repo root, in a Python >=3.10 environment:
pip install -e ".[dev]"        # installs vflank + ruff, mypy, pytest (needs the hatchling backend)
```

If you cannot install the build backend, you can still develop and test without
installing — `pyproject.toml` sets `pythonpath = ["src"]`, so:

```bash
python -m pytest               # tests resolve `vflank` from src/ automatically
PYTHONPATH=src python -m vflank.cli.app --help
```

Runtime dependencies: `typer`, `rich`, `pysam`, `pandas`. External data (not
bundled): an indexed reference FASTA (`.fai`), and optionally a directory of
gnomAD per-chromosome VCFs (`.bgz` + `.tbi`).

---

## 2. Project layout and the dependency rule

```
src/vflank/
├── core/     pure domain logic — NO file/CLI I/O
│   ├── chrom.py       chromosome notation detect/normalise
│   ├── variant.py     the Variant value object + validation
│   ├── flanks.py      FlankSource protocol, ReferenceFlankSource, mask_sequence
│   ├── popfreq.py     gnomAD resolve + parse_common_snp_positions + GnomadStore
│   ├── popfreq_api.py gnomAD GraphQL API source (GnomadApiSource)
│   ├── reference_api.py reference sequence via the UCSC API (ReferenceApiSource)
│   ├── consensus.py   BAM patient consensus (modes C/D) + insertion flagging
│   ├── fusion.py      Breakpoint/Fusion model + reverse-complement junction builder
│   └── skips.py       categorised skip-reason helpers
├── io/       file access
│   ├── maf.py         load/validate MAF (path or buffer), row -> Variant
│   ├── reference.py   ReferenceFasta + genome-build guard
│   ├── fasta.py       header sanitising + record formatting/writing
│   ├── breakpoints.py SV/fusion breakpoint-TSV reader (path or buffer)
│   ├── emit_primer3.py Primer3 Boulder-IO writer
│   └── report.py      TSV run-report writer
├── sources.py  reference/gnomAD source factories from config (validate + build)
├── pipeline.py the use case: iter_small/iter_fusion + collect + run_small/run_fusion
├── cli/      Typer commands (presentation only)
│   ├── app.py         root app, global -v/-q/--debug, version
│   ├── small.py       small-variant commands: run / inspect / list-vcf
│   ├── fusion.py      fusion command: run
│   ├── _bam.py        --bam/--bam-map resolver + ConsensusPolicy builder
│   └── _ui.py         parameter-echo panel
├── logging.py  shared Rich console + logger
└── errors.py   VflankError hierarchy
```

**Import direction is one-way:** `cli → pipeline → {sources, io} → core`. `core`
imports nothing from the layers above it; `pipeline` coordinates `sources`/`io`/
`core` but imports no `cli`/Rich/Typer. Keep it that way — it is what makes
`core` and `pipeline` unit-testable without pysam, pandas, or Typer.

---

## 3. Run the CLI

```bash
# Main pipeline: extract + mask flanks for every variant in a MAF
vflank small run variants.maf \
    --ref-genome /path/GRCh37.fasta \
    --pop-vcf-dir /path/gnomad_v2.1.1/ \  # local VCFs (--pop-source vcf, default)
    --genome-build hg19 \            # default; gnomAD v2.1.1. Use -g hg38 for GRCh38/v4.
    --pop-data genome \              # genome (default) | exome | both (union)
    --flank 200 \
    --output flanking_sequences.fasta \
    --report run_report.tsv          # optional: per-variant TSV + stats + skip breakdown

# No-download masking via the gnomAD API instead of local VCFs:
vflank small run variants.maf -r /path/GRCh37.fasta -g hg19 --pop-source api

# Fully no-download: reference from the UCSC API too (no local FASTA needed):
vflank small run variants.maf -g hg19 --ref-source api --pop-source api

# Preview MAF columns before a run (catches column-name mismatches)
vflank small inspect variants.maf

# Verify which gnomAD per-chromosome VCFs are found in a directory
vflank small list-vcf /path/gnomad_v4/ --genome-build hg38

# Global verbosity (before the subcommand): -v debug, -q quiet, --debug + tracebacks
vflank -v small run ...
```

`vflank small run --help` documents every option. Status output goes to
**stderr**; the FASTA goes to the `--output` file.

---

## 4. The quality gate

Run all three before considering any change done (this is enforced expectation,
see [CLAUDE.md](https://github.com/rhshah/vFlank/blob/main/CLAUDE.md)):

```bash
python -m ruff check src tests          # lint + import order + pyupgrade
python -m mypy src/vflank/core src/vflank/io src/vflank/pipeline.py
python -m pytest                        # unit + integration
```

Tests are split into `tests/unit/` (pure functions — no pysam needed) and
`tests/integration/` (build a tiny indexed FASTA and drive the real CLI via
Typer's `CliRunner`). Integration tests `pytest.importorskip("pysam")` so the
unit suite still runs in a minimal environment.

---

## 5. Use vflank as a library

### The high-level entrypoints (what a web service / notebook calls)

`vflank.pipeline.run_small` / `run_fusion` run the whole pipeline — build the
reference + gnomAD sources from options, load the input (a path **or** an open
buffer), orchestrate, close — and return a `RunResult` (records, per-variant
rows, categorised skips, counts, API-request tallies). No printing, no files
written; the caller decides what to render or persist.

```python
import io
from vflank.pipeline import run_small

maf = io.StringIO(open("variants.maf").read())     # a path works too
result = run_small(
    maf, genome_build="hg19",
    ref_source="api",                               # UCSC — no local FASTA
    pop_source="api",                               # gnomAD — no download
    flank=200, emit_primer3=True,
)
print(result.n_processed, result.n_skipped, len(result.records))
open("out.fasta", "w").writelines(result.records)
```

For incremental progress or custom accumulation, drop to the streaming
primitives `iter_small(df, ...)` / `iter_fusion(df, ...)` (generators of
`Processed | Skipped | Duplicate`) plus `collect(...)` — this is exactly how the
CLI drives its progress bar. See `docs/research/orchestration-extraction.md`.

### Composing the lower-level pieces

The building blocks are also usable directly. Minimal end-to-end example
against an indexed FASTA (no gnomAD; masked == raw):

```python
from vflank.core.variant import Variant
from vflank.core.flanks import ReferenceFlankSource
from vflank.io.reference import ReferenceFasta
from vflank.io import fasta as fasta_io

reference = ReferenceFasta("/path/GRCh38.fasta")          # requires .fai alongside
reference.check_build("hg38")                              # returns a warning string on mismatch, else None

source = ReferenceFlankSource(reference, gnomad=None, flank=200)

variant = Variant(chrom="7", start=140753336, end=140753336, ref="A", alt="T",
                  gene="BRAF", sample="P-001")
result = source.fetch(variant)                             # -> FlankResult
print(len(result.left), len(result.right), result.n_masked)

records = fasta_io.format_records(variant, result.upper(), variant.ref, variant.alt)
fasta_io.write_fasta("out.fasta", records)
reference.close()
```

Add gnomAD masking by passing a store:

```python
from vflank.core.popfreq import GnomadStore
gnomad = GnomadStore("/path/gnomad_v4/", "hg38")
source = ReferenceFlankSource(reference, gnomad=gnomad, flank=200, af_threshold=0.001)
# ... fetch variants ...
gnomad.close()
```

Drive from a MAF instead of constructing `Variant`s by hand:

```python
from vflank.io.maf import load_maf, parse_variant_row, MafColumns
cols = MafColumns()                                        # defaults to TCGA/MSK names
df = load_maf("variants.maf", cols)                        # validates required columns
for _, row in df.iterrows():
    variant, reason = parse_variant_row(row, cols)         # reason is a skip message or None
    if variant is None:
        continue
    result = source.fetch(variant)
```

The **pure kernels** need no files and are the easiest entry points to test or reuse:

```python
from vflank.core.chrom import normalise_chrom
from vflank.core.popfreq import parse_common_snp_positions
from vflank.core.flanks import mask_sequence

normalise_chrom("chr23")                                   # -> ("X", None)
parse_common_snp_positions(["chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.2"], 0.001)  # -> [100]
mask_sequence("ACGTACGT", 10, [11, 15])                    # region starts at 0-based 10 -> "NCGTNCGT"
```

---

## 6. Extending the codebase

### Add a new flank source (e.g. BAM consensus, modes C/D)

Implement the `FlankSource` protocol — anything with
`fetch(self, variant: Variant) -> FlankResult` works wherever
`ReferenceFlankSource` does:

```python
# src/vflank/core/consensus.py
from vflank.core.flanks import FlankResult
from vflank.core.variant import Variant

class ConsensusFlankSource:
    """Flanks from a sample BAM's consensus, masking het / low-coverage sites."""
    def __init__(self, reference, bam, *, flank=200, min_depth=20, het_frac=0.25):
        ...
    def fetch(self, variant: Variant) -> FlankResult:
        # build per-base consensus over the flank windows; fall back to
        # reference below min_depth; mask het/low-confidence positions to N.
        ...
```

Keep the coordinate math identical to `ReferenceFlankSource` (1-based MAF →
0-based half-open windows; variant interval excluded). Validate the consensus
against `samtools consensus` (see the `ddpcr-conventions` skill).

### Add a new population-frequency source

`get_positions(bare, start_0based, end_0based, af_threshold) -> list[int]` is the
masking interface `ReferenceFlankSource` depends on (duck-typed). Two
implementations already exist and are interchangeable:

- `core/popfreq.GnomadStore` — local gnomAD VCFs (`--pop-source vcf`). New VCF
  filename patterns go in `GNOMAD_PATTERNS` (keyed by build then `genome`/`exome`);
  reuse the pure `parse_common_snp_positions`.
- `core/popfreq_api.GnomadApiSource` — gnomAD GraphQL API (`--pop-source api`).
  Pure parser `parse_api_variants`; HTTP/throttle/clock are injected for offline
  tests. Honours the rate limit (region cache + throttle + backoff).

Both honour `--pop-data` via `kinds_for(...)`. To add another source (dbSNP,
1000G), implement the same `get_positions` and add a `--pop-source` branch.

### Add a new reference source

`ReferenceFlankSource` depends (duck-typed) on a reference exposing
`fetch(bare, start_0based, end_0based) -> str` (0-based half-open) plus
`check_build(declared) -> str | None`, `has_chr`, and `close()`. Two
implementations exist, selected by `--ref-source` via `cli/_reference.make_reference_source`:

- `io/reference.ReferenceFasta` — local indexed FASTA (`--ref-source file`,
  default); fingerprints the build by chr1 length.
- `core/reference_api.ReferenceApiSource` — UCSC getData/sequence API
  (`--ref-source api`, no download). Pure kernels `build_url` /
  `parse_sequence_response`; HTTP/throttle/clock injected for offline tests;
  trusts the requested build (no local sequence to fingerprint); exposes
  `request_count` for monitoring. See [research/genome-api.md](research/genome-api.md).

A short window off a contig end must return a **short string** (not raise) —
truncation is reported by the CLI, matching the pysam contract.

### Add a new CLI subcommand or sub-app

`cli/app.py` mounts sub-apps with `app.add_typer(...)`. Add a command to the
existing `small` app, or create a new module (e.g. `cli/fusion.py`) exposing
`app = typer.Typer()` and mount it. Conventions: status/Rich output via the
shared `console` (stderr); raise `VflankError` for user errors and let the
command wrapper print them; never `print()` data to stdout except the intended
output.

---

## 7. Invariants to preserve

- **Coordinates:** MAF is 1-based fully-closed; pysam is 0-based half-open. The
  flank formulas live in `ReferenceFlankSource.fetch` and `mask_sequence` — match
  them exactly in any new source.
- **No silent failures:** surface every error path (raise/log/report). Truncated
  flanks, absent contigs, and build mismatches are already handled this way —
  follow the pattern.
- **`core` stays pure and I/O-free** so the hot kernels remain testable and
  later portable to Rust.
- **Output contract:** two FASTA records per variant (raw + `Masked__`), variant
  shown literally as `[REF/ALT]`, only single-base population SNPs masked.

---

## 8. Contributing workflow

The repo follows **git-flow**: `main` holds tagged releases (`vX.Y.Z`) only,
`develop` is the integration branch. Initialise once with `git flow init`
(already configured here: prefixes `feature/ release/ hotfix/`, version-tag
prefix `v`). The installed git-flow is nvie 0.4.1, whose subcommands are
`feature`, `release`, `hotfix`, `support` — there is **no `bugfix` subcommand**
(that is the AVH edition); use a `feature` branch for fixes.

1. Start a feature off `develop`: `git flow feature start <name>`
   (equivalently `git checkout -b feature/<name> develop`).
2. Make the change; add/adjust tests in the matching `tests/` subtree.
3. Run the full quality gate (section 4) until green.
4. Re-read your diff for duplication, dead code, and unused symbols.
5. Commit, ending the message with:
   `Co-Authored-By: Claude Code <noreply@anthropic.com>`
6. Finish the feature back into `develop`: `git flow feature finish <name>`
   (or open a PR targeting `develop`).
7. **Releases:** `git flow release start X.Y.Z`, bump the version, then
   `git flow release finish X.Y.Z` — this merges into both `main` and
   `develop` and tags `vX.Y.Z`. Push with `git push --tags origin main develop`,
   then `gh release create vX.Y.Z` to publish (that triggers PyPI + GHCR).
   Urgent production fixes use `git flow hotfix start X.Y.Z` off `main`.

## 9. Continuous integration

Three workflows in `.github/workflows/`, each scoped to the git-flow event that
needs it:

| Event | `ci.yml` (lint · type · test · docs-build) | `docs.yml` (mike deploy) | `release.yml` (publish) |
|-------|:--:|:--:|:--:|
| PR → `develop` / `main` | ✓ | — | — |
| push `develop` | ✓ | `dev` alias | — |
| push `main` | ✓ | — | — |
| `vX.Y.Z` tag | — | `X.Y.Z` + `latest` (default) | — |
| GitHub Release published | — | — | PyPI + GHCR |

- **`ci.yml`** runs on pushes to `main`/`develop` and PRs targeting them. The
  `docs` job runs a plain `mkdocs build` (no deploy): it exits non-zero on
  config/plugin/**dependency** errors, so a missing-plugin breakage is caught in
  the PR instead of at deploy time.
- **`docs.yml`** publishes versioned docs with `mike`: the rolling `dev` alias
  from `develop`, and `X.Y.Z` + `latest` from a `vX.Y.Z` tag. A `docs-deploy`
  concurrency group serialises the gh-pages pushes.
- **`release.yml`** publishes to PyPI (OIDC Trusted Publishing) and GHCR when a
  GitHub Release is published. The release sequence (step 7) pushes the tag —
  which deploys versioned docs — then creates the Release — which publishes.
