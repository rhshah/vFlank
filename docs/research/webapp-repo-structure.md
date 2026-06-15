# Web app: in-repo vs. a separate `vFlank-webapp` — decision & plan

Status: **proposed** (2026-06-15). Decision input for: do we host the web tool
*inside* this repo, or stand up a separate **`vFlank-webapp`** repo? Builds on
[web-app-and-hosting.md](web-app-and-hosting.md) (product + hosting) and the
now-shipped `--ref-source api` (0.4.0), which removes the v1 reference
prerequisite that note flagged.

## Recommendation (TL;DR)

**A separate `vFlank-webapp` repo that depends on `vflank` from PyPI — but only
after one library-side change that belongs in *this* repo regardless: extracting
the per-variant orchestration out of the CLI into a reusable API.**

The web app should be a **downstream consumer** of vflank, exactly as Olivar /
Primer3 are downstream of the emit formats. Keep the library lean and
scope-bounded; keep web/JS/WASM concerns out of it.

The gating work is not the repo — it's the **orchestration extraction** below.
Decide the repo, but schedule the extraction first.

## The fact that decides it

The web-app note assumes v1 can "reuse the existing per-row loop as a library,
no new code." That is **not true today.** Verified:

- `cli/small.py` has **~110 Typer/Rich/console touch-points**; the per-variant
  loop, skip aggregation, report assembly, and summary all live inside
  `cli/small.py:_run` (and `cli/fusion.py:_run`), interleaved with progress bars
  and console output.
- There is **no presentation-free `pipeline`/`run` entrypoint** in `vflank`. A
  consumer today can import the *pieces* (`load_maf`, `parse_variant_row`,
  `ReferenceFlankSource`, `ReferenceApiSource`, `GnomadApiSource`,
  `format_records`) but would have to **re-implement the loop, skip handling, and
  dedup** — i.e. duplicate the core of `_run`.

So whichever repo the web app lives in, the library needs a reusable
orchestration function first. That requirement, not aesthetics, is the spine of
this plan.

## Step 0 — library prerequisite (lives in `vFlank`, any repo decision)

Extract a **presentation-free orchestration API** and leave the CLI a thin
presenter over it:

```python
# src/vflank/pipeline.py  (new; pure of Typer/Rich)
@dataclass
class RunResult:
    records: list[str]                 # FASTA records
    rows: list[dict]                   # per-variant detail (the report table)
    skips: dict[str, int]              # categorised skip counts
    stats: dict[str, object]           # provenance + outcome counters
    primer3: list[Primer3Record]       # when emit requested

def run_small(maf, *, reference, gnomad, flank, af_threshold, dedup,
              sample_filter=None, consensus=None, emit_primer3=False,
              progress=None) -> RunResult: ...
def run_fusion(sv_table, *, reference, gnomad, flank, ...) -> RunResult: ...
```

- `cli/small.py:_run` becomes: build sources → `run_small(...)` → render the
  table / write files. Rich progress is injected via an optional `progress`
  callback (default `None` = silent), so the library stays presentation-free and
  honours the `cli → io → core` dependency rule.
- Make `io.load_maf` accept a **path or file-like buffer** (one-line change) so a
  service needn't round-trip through a temp file.
- Publish a small, documented **public surface** (`vflank.pipeline.run_small/
  run_fusion`, the API sources, the emit writers). Add an `__all__` / a "Using
  vflank as a library" doc section (the DEVELOPER guide already starts this).

Wins beyond the web app: the orchestration becomes unit-testable **without**
`CliRunner`, and any future surface (Nextflow wrapper, notebook, the web app)
calls one function. This is the single highest-leverage refactor on the roadmap.

## Why a separate repo (the decision)

Once Step 0 lands, the web app is a thin consumer, and separation is the cleaner
architecture:

1. **Dependency hygiene.** The library's runtime deps are `typer/rich/pysam/
   pandas`. A web app adds `fastapi/uvicorn` (v1) and a **JS/WASM toolchain**
   (v2). Keeping those out of `vflank`'s `pyproject` preserves the "no
   speculative deps" discipline (CLAUDE.md) and keeps `pip install vflank` lean.
2. **Deploy lifecycle.** The library releases to **PyPI on tags**; the web app
   **auto-deploys to Render on push**. Different cadences, different triggers —
   entangling them means every webapp tweak risks the library's release flow and
   vice-versa.
3. **Frontend stack.** HTML/CSS/JS, a build step, and (v2) WASM assets +
   `node_modules` do not belong in a Python `src/`-layout package. v2 is a real
   frontend codebase, not a script.
4. **Scope & framing.** vflank is "the variant-aware masked-flank front-end";
   the web app is a *surface* that consumes it. A separate repo makes the
   patient-adjacent tool's (eventual, v2) PHI/compliance story self-contained.
5. **Contributors / issues.** Library users (pip/CLI, bioinformaticians) vs. web
   users are different audiences with different bug reports.

This mirrors how vflank already delegates **out-of-process** to Olivar/Primer3:
the boundary is a clean dependency edge, not a folder.

## When in-repo would win instead (honest counter)

- **If v1 stays a tiny, permanent single artifact** and you value one CI / one
  release / one issue tracker over separation — a `webapp/` subdirectory with its
  **own** `pyproject` and Render `rootDir: webapp` is a viable middle ground
  (monorepo, independent deploy). Costs: CI path-filtering, mixed Python-package
  + web concerns in one tree, and v2's JS/WASM still bloats it later.
- **If you expect to refactor library and service in lockstep constantly**,
  co-location reduces cross-repo PR friction — but Step 0 is designed precisely
  to make the contract stable, which removes most of that friction.

Net: the middle ground is defensible for v1-only; it gets worse exactly when v2
(WASM) arrives. Standing up the separate repo now avoids a later migration.

## `vFlank-webapp` plan (v1 = modes A/B, no PHI)

- **Stack:** FastAPI + uvicorn; a minimal server-rendered form (HTML + a little
  htmx/vanilla JS — **no build step** for v1). `pip install "vflank>=0.4"`
  (compatible-range pin).
- **Endpoints:** `GET /` (the 3-control form) and `POST /run` (multipart: mode,
  build, file) → JSON `{records, rows, skips, stats}`. Fixed sources:
  `--ref-source api` (UCSC, server-side — CORS is moot, see below) + gnomAD API.
- **Reuse:** call `vflank.pipeline.run_small/run_fusion`; **inherit all input
  validation** (the note's central point — now real, thanks to Step 0).
- **Service-only policy (not in the library):** a **≤10-record cap** post-parse
  (protect the shared instance / API rate limits); map `VflankError`/`MafError`
  → 4xx; surface per-row skips in the body. Deliberately not a cap in `run` (it
  would regress legitimate local batch users).
- **Hosting:** Render free tier (`render.yaml`), stateless, single instance — the
  natural home for the UCSC ~1 req/s throttle + an in-memory window cache.
  Cold-start after idle is acceptable for an internal/demo tool.
- **CI:** lint + test on PR; Render auto-deploys `main`. Its own (continuous)
  cadence, decoupled from PyPI.
- **Versioning:** the webapp pins a compatible `vflank` range and bumps when it
  needs a new library feature. The library never depends on the webapp.

## Phasing

| Phase | Scope | Repo pressure |
|---|---|---|
| **0** | Orchestration extraction + buffer `load_maf` + public API (in `vFlank`) | — (prereq) |
| **v1** | FastAPI service, modes A/B, UCSC+gnomAD APIs, no PHI | separate repo viable; monorepo-subdir tolerable |
| **v2** | BAM via client-side WASM (biowasm/Aioli) + JS/WASM port of the pure overlay; static Pages or hybrid | **decisively** separate (frontend codebase) |

The v2 pure-kernel port (`mask_sequence` + flank math + overlay) can be its own
small package the webapp vendors — another reason the boundary wants to be a
dependency edge, not a directory.

## Risks / open items

- **Step 0 sizing.** The extraction is the real effort; size it before
  committing dates. It is contained (move logic out of two `_run`s; the pieces
  are already pure) but touches the two busiest CLI files — do it behind the
  existing test suite (the `CliRunner` integration tests pin current behaviour).
- **UCSC CORS is moot for server v1.** Server-side calls aren't subject to
  browser CORS — already proven by the shipped `--ref-source api` working
  server-side. Only a static/Pages **v2** would need Ensembl (per the web-app
  note's correction table).
- **Two repos = two CI/release flows.** Accepted cost; mitigated by the stable
  public API from Step 0.
- **Render cold start / shared-IP throttle** — covered in the web-app note.

## Decision checklist

1. Approve **Step 0** (orchestration extraction) — needed regardless; do it next,
   in `vFlank`, as its own feature → a 0.5.0.
2. Choose repo: **separate `vFlank-webapp`** (recommended) vs. `webapp/` subdir
   (only if v1-forever).
3. After Step 0 ships, scaffold the chosen home with the v1 plan above.
