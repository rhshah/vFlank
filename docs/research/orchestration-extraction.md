# Orchestration extraction — methods & implementation plan

Status: **proposed** (2026-06-15). Prerequisite ("Step 0") from
[webapp-repo-structure.md](webapp-repo-structure.md): pull the per-variant
orchestration out of the CLI into a reusable, presentation-free `vflank` API so
any surface (web service, Nextflow, notebook) calls one function instead of
re-implementing the loop. This note covers the established methods and a
PR-sized plan grounded in the current code.

## The problem, concretely

`cli/small.py:_run` is ~380 lines that fuse six concerns:

| Concern | Examples in `_run` today |
|---|---|
| Presentation | `console.rule`, `echo_parameters`, ~30 `console.print`, `Progress`, the summary `Table` |
| Validation | genome-build, `validate_ref_source`, `validate_pop_options`, dir checks |
| I/O | `load_maf`, `write_fasta`, `write_report`, `write_primer3` |
| Source construction | `make_reference_source`, `make_pop_source`, the per-sample consensus cache + `_source_for` |
| **Orchestration** | the `for row in df` loop: parse → dedup → fetch → truncation → mask → format → primer3 → row/stat accumulation |
| Cleanup | `reference.close()`, `gnomad.close()`, closing cached consensus BAMs |

Only the **orchestration** row is the reusable core; everything else is shell.
`fusion.py:_run` has the same shape. There is no presentation-free entrypoint, so
a consumer must duplicate the loop, dedup, skip-aggregation, and stats.

## Current best methods (the "how")

This is a textbook **CLI-to-library** extraction. The relevant, current patterns:

1. **Functional core, imperative shell** (Bernhardt). Push decision-making into a
   core that takes values and returns values; keep I/O and rendering in a thin
   shell. Here: the loop becomes the core; `console`/`Progress`/file-writes stay
   in the shell.
2. **Hexagonal / ports & adapters.** vflank is already **half-way**: `FlankSource`,
   the duck-typed `get_positions`, `ReferenceFasta`/`ReferenceApiSource`,
   `GnomadStore`/`GnomadApiSource` are ports with swappable adapters. The missing
   piece is the **use-case layer** (the orchestration) that sits above the ports;
   extracting it *completes* the hexagon. CLI and web become two adapters over
   one use case.
3. **Command–Query Separation / return-don't-print.** The core returns a result
   object; it never prints or writes files. Callers decide what to render/persist.
4. **Effects injected, not imported** (the "Sans-IO" discipline). The core
   receives its sources (reference, gnomAD, BAM resolver) as parameters; it does
   not construct them or reach for the network/filesystem itself. Construction is
   the shell's job (config → adapter).
5. **Progress as a stream, not a print.** Three real options:
   - *Callback*: pass `on_progress(done, total)`. Works, but the core still
     "knows about" progress and you thread a parameter through.
   - *Generator / event stream* (**recommended**): the orchestration is an
     iterable that **yields one outcome per variant**. The shell drives progress
     by counting yields and has zero coupling the other way. Bonus: streaming —
     incremental file writes, bounded memory, and early results for a web
     response. It composes with Rich directly: `for o in Progress(...).track(
     iter_small(...), total=n)`. This is the idiomatic modern-Python choice
     (mirrors `rich.progress.track`, which wraps an iterable).
   - *Observer Protocol*: an object with `on_variant`/`on_skip`. More ceremony;
     rarely worth it over a generator in Python.
6. **Diagnostics via `logging`** (already in place: `log.warning` for missing
   BAMs). The core logs; the shell configures handlers. No `print` in the core.
7. **Characterization tests as the safety net.** Before moving code, treat the
   existing `CliRunner` integration tests (which assert exact FASTA + report
   bytes) as a **golden master**: don't touch them, keep them green through every
   step — that *is* the proof the refactor preserved behaviour.
8. **Incremental, behaviour-preserving moves** (strangler-style). Extract in
   small PRs that each keep the gate green, rather than one big rewrite.

## Target architecture

Introduce a **use-case layer** between the CLI and io/core:

```
cli/  (adapter: flags → config, render, write files, Rich progress)
  │
pipeline.py  (use case: the per-variant orchestration — NO Typer/Rich/print)
  │
io/   (load_maf, writers)      core/  (sources, flanks, mask, fusion — pure)
```

`pipeline` may use `io` and `core` (it coordinates them); it must not import
`cli` or `rich`/`typer`. This extends the existing one-way rule to
`cli → pipeline → io → core`.

## The design

**A generator of discriminated outcomes + an eager collector.**

```python
# vflank/pipeline.py  (new)
@dataclass(frozen=True, slots=True)
class Processed:
    variant: Variant
    records: list[str]            # FASTA (raw + masked)
    detail: dict                  # the report row
    primer3: Primer3Record | None
    used_consensus: bool
    truncation: str | None        # warning text, or None
    flagged: bool

@dataclass(frozen=True, slots=True)
class Skipped:
    ref: str                      # "row 12 …" identity
    reason: str

@dataclass(frozen=True, slots=True)
class Duplicate:
    variant: Variant

Outcome = Processed | Skipped | Duplicate

def iter_small(
    df, *, cols, reference, gnomad, flank, af_threshold, dedup, uppercase,
    bam_resolver=None, policy=None, emit_primer3=False,
) -> Iterator[Outcome]:
    """Yield one Outcome per MAF row. Owns the per-sample consensus cache and
    closes it on exit (try/finally) — pure of Typer/Rich/file-writes."""
    cache: dict[str, object] = {}
    try:
        for row_idx, row in df.iterrows():
            ...  # the body lifted verbatim from _run, `yield` instead of append
    finally:
        for src in cache.values():
            if hasattr(src, "bam"): src.bam.close()

@dataclass
class RunResult:
    records: list[str]
    rows: list[dict]
    skips: dict[str, int]
    skip_messages: list[str]
    truncations: list[str]
    primer3: list[Primer3Record]
    stats: dict[str, object]

def collect(outcomes: Iterable[Outcome]) -> RunResult:
    """Accumulate a stream into a RunResult (records, rows, counts, stats)."""

def run_small(config: SmallConfig) -> RunResult:
    """Batteries-included: validate → build sources → iter → collect → close.
    The entrypoint a web service / notebook calls."""
```

- **CLI** keeps its UX: `df = load_maf(...)`; build sources (adapters);
  `for o in track(iter_small(df, …), total=len(df)): accumulate`; then write
  files and render the summary `Table`. Progress falls out of iteration; nothing
  in `pipeline` knows Rich exists.
- **Web** calls `run_small(config)` and serialises `RunResult` to JSON.
- **Ownership rule for cleanup:** whoever constructs a source closes it. The
  caller builds & closes `reference`/`gnomad`; `iter_small` builds & closes the
  consensus cache (in `finally`). `run_small` builds everything, so it closes
  everything.

**Validation & source factories move to the library** so both adapters inherit
them (today they live in `cli/_reference.py`, `cli/_masking.py`, and inline
checks in `_run`):

```
vflank/config.py     SmallConfig/FusionConfig dataclasses + validate()
vflank/sources.py    make_reference_source, make_pop_source, build_consensus_resolver
```

These are already presentation-free (they raise `VflankError`, no Rich), so this
is a move, not a rewrite.

## What moves vs. stays

| Piece | New home |
|---|---|
| the `for row in df` loop, dedup, `_source_for`, consensus cache, truncation check, stat counters | `pipeline.iter_small` |
| `RunResult` assembly (records, rows, skip/dup counts, stats) | `pipeline.collect` |
| `validate_*`, genome-build guard, dir checks | `config.validate` |
| `make_reference_source`, `make_pop_source`, consensus resolver | `sources.py` |
| `load_maf`, `write_fasta`, `write_report`, `write_primer3` | stay in `io`; **called by the shell** (CLI) or by `run_small` |
| `echo_parameters`, every `console.print`, `Progress`, summary `Table` | **stay in `cli`** |

## Implementation plan (PR-sized, gate green at each step)

1. ✅ **PR1 — `pipeline.iter_small` + `collect` (small path), behaviour-preserving.**
   Lifted the loop and `_source_for`/cache out of `_run` into `iter_small`; `_run`
   consumes the generator and keeps *all* presentation, writes, and cleanup of
   reference/gnomAD. No output changes. Direct pysam-free unit tests for
   `iter_small`. `pipeline.py` added to the mypy gate. **Done.**
3. ✅ **PR3 — fusion.** `iter_fusion`, same treatment; `Processed` generalised
   (the unused `variant` field dropped) so `collect` is shared. **Done.**

**Deferred until the web service exists (no speculative code — CLAUDE.md):** the
following have **no caller** until `vFlank-webapp` is scaffolded, so building them
now would be unused code. They land *with* that work, or when the CLI is routed
through `run_small` to give it a real caller:

2. **PR2 — `config.py` + `sources.py` + `run_small`.** Move validators and the
   source factories out of `cli/`; add `SmallConfig`/`FusionConfig` +
   `run_small`/`run_fusion(config) -> RunResult` (build sources → iterate →
   collect → close → assemble stats + request counts). This is non-speculative
   the moment the CLI `_run` routes through it (eliminating the source-build +
   close duplication between CLI and web).
4. **PR4 — buffer input.** `load_maf` / `load_sv_table` accept a path *or* a
   file-like buffer so a service needn't touch a temp file.
5. **PR5 — public API + docs.** `__all__`, a "Using vflank as a library" section
   for `run_small`/`run_fusion` + `RunResult`. Ship as **0.5.0**.

The keystone (PR1+PR3) is done and `develop` is shippable; the rest is pulled by
the consumer, not pushed ahead of it.

## Testing strategy

- **Before:** the existing `tests/integration/*` (`CliRunner`, exact FASTA +
  report asserts) are frozen as the golden master — unchanged, green throughout.
- **After PR1+:** add `tests/unit/test_pipeline.py` exercising `iter_small`/
  `collect`/`run_small` directly with a tiny indexed FASTA — faster and more
  granular than driving the CLI, and the path a web test will reuse.
- **Parity check:** one test asserts `run_small(config).records` equals the FASTA
  the CLI writes for the same inputs (the two adapters agree).

## Risks & edge cases to preserve exactly

- **Dedup key** includes the sample in BAM mode (one record per variant×sample).
- **Per-row fetch errors** become a Skip with a message — never abort the run.
- **Truncation semantics:** `exp_left = min(flank, start-1)`; short right (or
  shorter-than-allowed left) flags truncation but still emits the record.
- **`uppercase`** applies to the flanks *and* ref/alt together.
- **Cleanup order / ownership** as above — verify no leaked BAM handles via the
  `finally` path (test early-exit).
- **`request_count`** for the reference/gnomAD APIs is read off the sources after
  iteration and folded into `stats` — keep surfacing it in the report.
- **Sample filtering** stays a DataFrame pre-filter in a shared helper (it needs
  the sample column and is cheap), applied before `iter_small`.

## Effort

Moderate, risk concentrated in PR1/PR3 (the two busiest files). The pieces being
moved are already pure functions and protocols, and the golden-master tests pin
behaviour — so this is disciplined cut-and-lift, not redesign. The payoff is
outsized: orchestration becomes unit-testable without `CliRunner`, and it is the
single unlock for the web service, a Nextflow wrapper, and notebook use.
