# Single-variant mode, web UI, and hosting — design note

Exploratory design note for turning vflank into an interactive **single-variant**
web tool, and where to host it. Status: **proposed** (nothing built yet). This
note exists to make the trade-offs reviewable before committing, in the spirit
of `gnomad-api.md`. Verified facts are dated; open questions are flagged.

## Summary / recommendation

The three ideas — a single-variant mode, a UI, and WebAssembly for BAM — stack
into one coherent product: **a single-variant web tool where the reference base
and population frequencies come from public APIs, and the only large/private
input (the patient BAM) is computed client-side in the browser and never
uploaded.** Batch/cohort work stays where it is today: server-side CLI + local
files, fanned out by Nextflow.

Recommended path, in order:

1. **Single-variant mode** (`small one` / `fusion one`) — the request/response
   primitive every other piece needs. Low risk: the core is already per-variant.
2. **API-backed sources** — gnomAD API exists (`--pop-source api`); add a
   reference API source (see `genome-api.md`, prefer UCSC). Together these make a
   stateless server possible (no local FASTA/VCF).
3. **UI** — a thin layer over the single-variant endpoint.
4. **WASM-BAM** — `samtools` via biowasm/Aioli, client-side, for the
   patient-consensus path; keeps PHI in the browser.

**Hosting recommendation:** for the fastest path that reuses the existing Python
core, host a small FastAPI service on **Render** (or Fly.io / Railway / Cloud
Run) with the BAM step running client-side in WASM. A fully static **GitHub
Pages** deployment is achievable *only* if the Python kernel is ported to
JS/WASM and every external API allows browser CORS — a larger lift, best as an
end-state or a public demo. A **hybrid** (Pages frontend + a small API backend)
is the pragmatic middle.

## 1. Single-variant mode — the primitive

Both CLI entry points are file-first today, but the core is already
single-variant:

- `cli/small.py` loop: `maf_io.parse_variant_row(row, cols)` → one `Variant` →
  `source.fetch(variant)` → `FlankResult`.
- `cli/fusion.py` loop: `bp_io.parse_fusion_row(row, cols)` → one `Fusion` →
  `core.fusion.build_junction(...)`.

A "one variant" entry constructs a single `Variant`/`Fusion` from arguments
instead of a parsed row and calls the **identical** path:

```
vflank small one  --at 17:7579472 --ref C --alt G   -r REF -g hg19 [--pop-source api]
vflank fusion one --bp1 9:133589268:+ --bp2 22:23632600:+  -r REF -g hg19
```

Design constraints:

- **No duplicated orchestration.** Factor the shared "build source → fetch →
  format → (report) one record" out of `_run` so file-mode and one-mode share
  it. This is required by the repo's no-duplication discipline, not optional.
- This is the exact shape an HTTP endpoint and a UI form need: one variant in,
  raw + `Masked__` out, **no file upload**.
- Scope guardrail (per CLAUDE.md): it emits the masked target sequence and the
  Olivar/Primer3 emit formats; it does not become a primer designer.

## 2. API-backed sources — what makes a stateless server possible

- **gnomAD**: already implemented (`core/popfreq_api.GnomadApiSource`,
  `--pop-source api`). Rate limit (~10 req/60 s) rules out bulk but is a
  non-issue for interactive single-variant use (~1 request per variant window).
- **Reference**: proposed `ReferenceApiSource` mirroring `ReferenceFasta.fetch()`
  — see `genome-api.md`. UCSC preferred (0-based half-open coords + `hg19`/`hg38`
  genome names match the codebase exactly; no coordinate translation).

With both, a single-variant server needs **no local FASTA/VCF** — it can be a
small stateless container. Reference sequence is immutable, so the
reproducibility caveat that keeps the gnomAD API *optional* does not apply to the
reference API.

## 3. WebAssembly for BAM — the one input that can't be an API

Reference and gnomAD can be remote calls; **patient reads cannot** — BAMs are
large and are PHI. You can't (and shouldn't) upload a multi-GB BAM to a hosted
service. WASM resolves this tension: run `samtools` **client-side in the
browser**, so reads never leave the user's machine.

**Tooling:** [biowasm](https://biowasm.com/) compiles `samtools` (and htslib,
`bcftools`, etc.) to WebAssembly; [Aioli](https://github.com/biowasm/aioli) runs
them in a WebWorker with a virtual filesystem and **lazy byte-range reads of
local files**.

Why this fits vflank specifically:

- **Tiny windows.** vflank reads only ±`flank` bp (default 200) around the
  variant. With the `.bai` index, samtools does a *ranged* read — a few KB of
  BGZF blocks, not the whole BAM. Aioli's lazy local-file reads are the best case
  for this access pattern; the whole BAM is never loaded into the tab.
- **Engine already pluggable.** `bam-consensus.md` defines
  `--bam-engine {samtools,pileup}` and delegates base-calling to
  `samtools consensus`. A browser/WASM path is the *same samtools engine on a
  different host*, not a new algorithm.
- **Overlay already pure.** `apply_lowcov_overlay` / `mask_sequence` /
  `normalise_chrom` were deliberately kept pure (CLAUDE.md: portable to Rust).
  Porting the small overlay to JS or Rust→WASM is the move the design already
  anticipated. Split: `samtools-wasm` does base-calling; the ported pure overlay
  does the gnomAD-low-coverage layering and insertion masking.

Caveats / things to confirm:

- **Confirm the WASM samtools build exposes the `consensus` subcommand.** It
  needs samtools ≥ 1.16 (where `consensus` was introduced); biowasm ships 1.17/
  1.18 (verified below), so the version is sufficient — confirm the compiled
  module includes the subcommand at runtime.
- You'd run `samtools-wasm` + a ported overlay in the browser, **not** vflank's
  Python. The alternative — all of vflank via **Pyodide** (CPython-on-WASM) —
  founders on pysam (htslib C wrapper, not readily available in Pyodide). Don't
  lead with Pyodide; the biowasm-engine + ported-overlay route is lighter and
  more proven.

## 4. Hosting options

Two fundamentally different shapes, driven by **where the kernel runs**:

| Host | Model | Runs the Python core? | Cost / ops | Fit |
|---|---|---|---|---|
| **GitHub Pages / Netlify / Cloudflare Pages** | Fully static | **No** — kernel must be ported to JS/WASM (or Pyodide) | Free, zero infra, no server to scale | Public demo / truly serverless **end-state**; biggest lift |
| **Render / Fly.io / Railway / Cloud Run** | Small server | **Yes** — keeps existing FastAPI + Python core | Cheap; free tiers spin down (cold starts) | **Fastest path**; reuses all current code |
| **Pages frontend + small API backend** | Hybrid | Yes (backend only) | Free frontend + cheap backend | Pragmatic middle for a public-facing tool |

Decision drivers:

- **Static-only (Pages) requires browser CORS on every external call.** Verified:
  the gnomAD API sends `Access-Control-Allow-Origin: *`, so direct browser calls
  work. **Open:** UCSC/Ensembl reference-API CORS (below). If a reference API
  blocks cross-origin, a static site needs a tiny proxy — which means a server,
  defeating Pages-only.
- **Static-only also requires porting the kernel** (flank math + masking) to
  JS/WASM. That's real work but keeps the pure-function discipline honest.
- **A server (Render et al.) sidesteps both:** keep the Python kernel; server-side
  API calls avoid CORS entirely; the server stays stateless because reference +
  gnomAD are APIs.
- **The BAM step is client-side WASM regardless of host** — so PHI never touches
  the chosen platform. That's a strong property when hosting patient-adjacent
  tooling on a third-party PaaS: no patient reads on Render/Pages, ever.

## Verified facts (2026-06-14) and open questions

Verified:

- **gnomAD GraphQL API allows cross-origin browser requests.** Its server applies
  `app.use(cors())` (default config → `Access-Control-Allow-Origin: *`) in
  `broadinstitute/gnomad-browser` `graphql-api/src/app.ts`. → A static browser
  app can call gnomAD directly.
- **biowasm ships samtools 1.17 and 1.18** (both ≥ 1.16), so the `consensus`
  subcommand is available in principle.

Open (verify before committing to a static/Pages build):

- **UCSC / Ensembl reference-API CORS** — does the chosen reference API send
  permissive CORS headers for browser calls? If not, a reference proxy (small
  server) is required, and Pages-only is off the table.
- **Runtime confirmation** that the biowasm samtools module actually exposes
  `consensus` (not just that the version supports it).
- **Local BAM range-read ergonomics** in Aioli for a user-picked `.bam` + `.bai`
  (File System Access API vs. file slicing) on the target browsers.

## Scope guardrail

Whatever the host, vflank stays the variant-aware masked-flank front-end: it
serves target sequences and Olivar/Primer3 emit formats. A web surface does not
justify growing a primer/probe design algorithm into this package (see
ARCHITECTURE.md scope boundary).

## Sources

- gnomAD CORS: [`graphql-api`](https://github.com/broadinstitute/gnomad-browser/tree/main/graphql-api),
  [gnomAD API help](https://gnomad.broadinstitute.org/help/how-do-i-query-a-batch-of-variants-do-you-have-an-api)
- WASM tooling: [biowasm](https://github.com/biowasm/biowasm),
  [Aioli](https://github.com/biowasm/aioli)
</content>
</invoke>
