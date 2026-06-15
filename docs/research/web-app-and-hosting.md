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

Recommended path, in order (note the dependency: the service can't be fileless
on Render without the API sources, so they come **first**):

1. **API-backed sources — the v1 prerequisite.** Render has no local FASTA/VCF
   (hg19/hg38 are ~3 GB each; the disk is ephemeral), so the reference and
   population data *must* come from APIs. gnomAD already exists
   (`--pop-source api`); the one new core piece v1 needs is a **reference API
   source** (`ReferenceApiSource`, see `genome-api.md`). **Reference choice is
   host-dependent** (verified below): **UCSC** behind a server (nicer coords, no
   CORS needed — our case), **Ensembl** for a static/browser app (CORS-safe).
2. **Thin web service over the existing pipeline** — reuse `io/maf.load_maf` +
   the existing per-row loop as a *library* (no new CLI command, no refactor),
   run with `--ref-source api --pop-source api` so it touches no local files.
   The service owns a small-input cap and HTTP error mapping; **all existing
   input validation is inherited**. See §1.
3. **UI** — a thin 3-control front-end over the service (see §1, "The v1 form").
4. **WASM-BAM** — `samtools` via biowasm/Aioli, client-side, for the
   patient-consensus path; keeps PHI in the browser (v2).

A single-variant CLI mode (`small one` / `fusion one`) and the `make_variant`
refactor are **not required for v1** — they become worthwhile only if a
structured form UX is later preferred over file input (see §1, "Deferred").

### v1 scope decision — ship without BAM (steps 1–3)

**v1 = steps 1–3; BAM (step 4) is deferred to v2.** This is the cleanest first
release because the BAM is the *only* thing that forces the hard parts — WASM /
Aioli / the ported overlay, client-side file handling, the BAM-input spectrum,
**and the entire PHI/compliance story**. Dropping it for v1 means:

- v1 = single variant (SNV or SV) + reference API + gnomAD API → masked flanks,
  i.e. exactly the existing **modes A/B** (the CLI's default with no `--bam`). No
  new science.
- **No patient data ever touches the app** — only public reference + gnomAD data
  flow through it, so hosting carries no PHI liability anywhere.
- A fully static **GitHub Pages** build is genuinely viable for v1: gnomAD +
  Ensembl are both CORS-safe (verified), so the only remaining lift is porting
  the small pure kernel (`mask_sequence` + flank math) to JS — no WASM blocker.
- **BAM (modes C/D) becomes a self-contained v2** that adds the patient-consensus
  path and picks from the BAM-input options (below) then, without gating v1.

**Hosting recommendation:** for the fastest path that reuses the existing Python
core, host a small FastAPI service on **Render** (or Fly.io / Railway / Cloud
Run). With BAM deferred (v1), there is no PHI and a fully static **GitHub Pages**
build is also viable once the small pure kernel is ported to JS. When BAM lands
(v2), the client-side WASM path keeps patient reads in the browser; a **hybrid**
(Pages frontend + a small API backend) is the pragmatic middle.

### v1 hosting decision — Render free tier

**Decided: FastAPI on Render's free tier.** Chosen over the static/Pages route
because it reuses the existing Python core with **no JS port** (see §"the JS
port" reasoning), and the free tier is ample because v1 is not compute- or
IO-heavy: one request = two outbound API calls (reference + gnomAD) + a few
hundred bytes of string masking. No local FASTA/VCF, no pandas/pysam, no large
in-memory data, nothing persisted — a textbook fit for 512 MB / shared CPU and
an ephemeral filesystem.

Plan around these free-tier traits (none is a blocker for an internal/demo tool):

- **Spin-down → cold start.** The service sleeps after ~15 min idle; the first
  request after that takes tens of seconds to wake. Acceptable for v1; if it
  matters later, add a keep-warm ping or move to a cheap paid tier.
- **Latency is the external APIs, not Render.** Users wait on the UCSC + gnomAD
  round-trips. A small in-memory cache covers repeat lookups.
- **Single instance = one IP.** The free tier is one process, which is the
  natural place for the in-memory cache *and* the UCSC ~1 req/s throttle — one
  cache, one rate limiter, no coordination across instances.
- **Stateless, no disk** — matches the ephemeral filesystem; nothing to persist.

## 1. v1 implementation — reuse the existing pipeline (no new command)

**Decision:** v1 is a thin web service that calls the existing pipeline as a
**library** — no `small one` / `fusion one` CLI command, no `make_variant`
refactor, and **no new input validation**. The service is the only new code.

Why this is enough: the core is already per-variant, and the file path already
validates input at both levels, so the service inherits all of it:

- `io/maf.load_maf` — unreadable / empty file and **missing required columns** →
  `MafError` (with a remap hint); optional label columns defaulted.
- `io/maf.parse_variant_row` (per row) — `normalise_chrom`, numeric-position
  coercion, `validate_coordinates`, `validate_allele` → a skip reason the loop
  aggregates. Fusions: `io/breakpoints.parse_fusion_row` (chrom / pos / strand).

The validators (`validate_coordinates`, `validate_allele`, `normalise_chrom`)
are already standalone pure functions; `parse_*_row` is just pandas extraction
around them. So **no input validation is written or duplicated** for v1.

What the service adds — at its own layer, not in core:

- **Small-input cap (UI / service only).** Reject inputs over ~10 records
  *after* parse, with a clear "≤10 records in the hosted UI" message. This is a
  hosting *policy* guard (protect the shared instance / the API rate limit),
  deliberately **not** placed in `run`: a hard cap there would regress
  legitimate local-file batch users (no rate limit on the local path; the
  existing soft warning at >50 in `small.py` is already tied to
  `--pop-source api`, the correct condition).
- **HTTP mapping.** Translate `MafError` / `VflankError` → 4xx; surface per-row
  skip reasons in the response body.
- **Input plumbing.** `load_maf` takes a path, so the service writes the
  uploaded text to a temp file (fine on Render's ephemeral FS) — or a one-line
  refactor lets `load_maf` accept a file-like buffer.

Scope guardrail (per CLAUDE.md): the service emits the masked target sequence and
the Olivar / Primer3 emit formats; it does not become a primer designer.

### The v1 form (3 controls)

The UI is a **parameter + upload** form — not a structured per-variant form — so
it maps one-to-one onto existing pipeline parameters and needs no `make_variant`:

| Form control | Maps to | Notes |
|---|---|---|
| **Small vs Fusion** | `small run` vs `fusion run` pipeline | both paths already exist |
| **GRCh37 vs GRCh38** | `--genome-build hg19` / `hg38` | both builds supported end-to-end |
| **Upload file** | MAF (small) or breakpoint TSV (fusion) | service writes to a temp file → `load_maf` / `load_sv_table` |

Everything else takes a default hidden in v1 (flank ±200, AF threshold 0.001,
pop-data genome) and is exposed later if needed. Sources are fixed for the
hosted service: `--ref-source api` (UCSC) + `--pop-source api`.

### Deferred (only if a structured-variant form UX is later preferred)

A browser **form** (type `chr/pos/ref/alt`, or two breakpoints) is a nicer
single-variant UX than uploading a tiny MAF, but it needs structured input that
doesn't flow through pandas. That — and only that — is when the following earn
their place, and even then they **reuse the same validators**, not new logic:

- `core/variant.make_variant(...)` / `core/fusion.make_fusion(...)` — thin, pure
  orchestration of the existing validators minus pandas; `parse_*_row` delegates
  to them (removing the pandas/validation coupling).
- Optional `small one` / `fusion one` CLI commands as terminal conveniences over
  the same builders.

This keeps v1 minimal while leaving a clean, non-duplicating path to the form UX
if/when it's wanted.

## 2. API-backed sources — what makes a stateless server possible

- **gnomAD**: already implemented (`core/popfreq_api.GnomadApiSource`,
  `--pop-source api`). Rate limit (~10 req/60 s) rules out bulk but is a
  non-issue for interactive single-variant use (~1 request per variant window).
- **Reference**: `ReferenceApiSource` mirroring `ReferenceFasta.fetch()` — see
  `genome-api.md`. UCSC preferred (0-based half-open coords + `hg19`/`hg38`
  genome names match the codebase exactly; no coordinate translation). **This is
  the one new core piece v1 requires** — without it the hosted service has no
  reference, since Render carries no local FASTA. Wire it behind a
  `--ref-source file|api` flag mirroring `--pop-source`.

With both, the hosted service needs **no local FASTA/VCF** — it can be a small
stateless container. Reference sequence is immutable, so the reproducibility
caveat that keeps the gnomAD API *optional* does not apply to the reference API.

**Build-guard consequence.** Today `ReferenceFasta` fingerprints the build by
chr1 length to catch hg19/hg38 mix-ups. With an API reference there is no local
FASTA to fingerprint, so the service instead **trusts the build the form
selects** (the API serves exactly that assembly). `ReferenceApiSource` should
still pass the selected build through to the API and surface a clear error if
the API rejects it — no silent fallback.

## 3. BAM input — deferred to v2 (the one input that can't be an API)

**Deferred to v2** (see the v1 scope decision above). Captured here so the
research isn't lost. Reference and gnomAD can be remote calls; **patient reads
cannot** — BAMs are large and are PHI. The right axis for choosing is *where the
reads live and who sees them*, with file size secondary. Single-variant scope
helps every option: you only ever need one ±flank window, never the whole BAM.

### BAM-input options (v2)

*Reads never leave the user's machine:*

1. **Client-side WASM (biowasm/Aioli)** — local file picker, `.bai` ranged read,
   consensus in the tab. Strongest privacy; works even on static Pages. Cost:
   WASM engine + the pure overlay ported to JS/Rust. (Detailed below.)
2. **Client-side slice → upload only the window (hybrid)** — browser extracts the
   ±flank mini-BAM (`samtools view -b region` in WASM, a few KB) and uploads only
   that slice to a backend running the **real pysam `BamConsensusSource`**.
   Reuses the Python core verbatim; KB not GB; only reads near one locus leave
   the machine.

*Reads (or a slice) transit your server:*

3. **Server-side ranged fetch from a user-provided URL** — pysam/htslib opens a
   BAM/CRAM over `https://`/`s3://` and reads only indexed byte ranges; the file
   stays in the user's bucket (signed URL). IGV-style. `BamConsensusSource` works
   with a URL instead of a path. Needs the `.bai` alongside + access.
4. **GA4GH htsget** — standardized authorized, region-scoped read streaming.
   Purpose-built, but needs an htsget server to exist (most labs don't run one).
5. **Full server-side upload** — simplest to build (reuses everything), worst on
   PHI and size; for single-variant, options 2–3 dominate it.

*Reads stay in the secure pipeline:*

6. **Pre-compute, serve the derived sequence** — run consensus offline
   (CLI/Nextflow) in the secure environment; the web tool consumes only the
   already-masked per-variant sequence. Zero BAM in the app; not interactive for
   new variants.

*Ruled out:* **Pyodide + pysam** in-browser (htslib C not readily available in
Pyodide) — hence WASM uses samtools-wasm + a ported overlay, not Python.

*Cross-cutting:* **CRAM** (smaller, same indexed ranged access) helps 2–4.
**Host shape constrains the set** — static Pages can do #1 (and #3 if the bucket
sends CORS); a server unlocks #2–#6.

### Option 1 in detail — WASM resolves the privacy tension

WASM lets you run `samtools` **client-side in the browser**, so reads never leave
the user's machine.

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

- **Static-only (Pages) requires browser CORS on every external call.** Verified
  (below): gnomAD and **Ensembl** both send `Access-Control-Allow-Origin: *`, so
  a static app can call both directly. **UCSC apparently does not** set CORS for
  browser use — so a Pages app must use Ensembl for the reference, or proxy UCSC
  (a proxy means a server, defeating Pages-only).
- **Static-only also requires porting the kernel** (flank math + masking) to
  JS/WASM. That's real work but keeps the pure-function discipline honest.
- **A server (Render et al.) sidesteps both:** keep the Python kernel; server-side
  API calls avoid CORS entirely; the server stays stateless because reference +
  gnomAD are APIs.
- **The BAM step is client-side WASM regardless of host** — so PHI never touches
  the chosen platform. That's a strong property when hosting patient-adjacent
  tooling on a third-party PaaS: no patient reads on Render/Pages, ever.

## Verification log (2026-06-14)

**Environment limitation, stated up front.** This session's network policy
**blocks direct egress** to the genome/annotation API hosts — `curl` to
`api.genome.ucsc.edu`, `rest.ensembl.org`, and `gnomad.broadinstitute.org` all
returned `HTTP 403` (an egress proxy block, not the APIs themselves);
`raw.githubusercontent.com` was reachable (`301`). So CORS could **not** be
confirmed by capturing live response headers here. The findings below are from
**authoritative source/docs**, which is conclusive for the source-backed items
and strong (but not a live header capture) for the doc-backed ones. The one
cheap final step — open a browser, run `fetch()` from a throwaway origin, read
the `Access-Control-Allow-Origin` header in devtools — is noted per item.

### CONFIRMED — gnomAD API is browser-CORS-safe

- **Evidence (source):** `broadinstitute/gnomad-browser`,
  `graphql-api/src/app.ts` imports `cors` and calls `app.use(cors())` with no
  arguments. The `cors` package's default config emits
  `Access-Control-Allow-Origin: *`.
- **Conclusion:** a static browser app can POST GraphQL to
  `https://gnomad.broadinstitute.org/api` cross-origin. Rate limit (~10/60 s)
  still applies but is irrelevant for single-variant interactive use.
- **Confidence:** high (server source). Live header check optional.

### CONFIRMED — Ensembl REST is browser-CORS-safe

- **Evidence (docs):** the official Ensembl wiki page *CORS And JSONP*
  (`Ensembl/ensembl-rest`) states Ensembl REST returns
  `Access-Control-Allow-Origin: *` whenever an `Origin` request header is sent,
  and calls CORS "the best way to access data in Ensembl REST from a browser."
- **Conclusion:** Ensembl REST is usable directly from a static browser app.
  Note: GRCh37 lives on `grch37.rest.ensembl.org` (separate host); 1-based
  coordinates need a `+1` conversion vs. our 0-based half-open interface.
- **Confidence:** high (vendor docs). Live header check optional.

### LIKELY NOT browser-CORS-safe — UCSC REST API

- **Evidence (community + docs):** a UCSC `genome` mailing-list thread ("Our web
  application got blocked by api.genome.ucsc.edu") and multiple community
  write-ups indicate browser apps hit CORS blocking and/or aggressive
  rate-limit blocking against `api.genome.ucsc.edu`, with a **server-side proxy**
  given as the standard workaround. The UCSC API help page itself does not
  advertise CORS support, and recommends ~1 request/second.
- **Conclusion / impact:** UCSC is **not** a safe direct-from-browser reference
  source. This **changes the reference-API choice by host** (see correction
  below). It remains fine behind a server (server-side calls need no CORS).
- **Confidence:** medium-high (community + absence of vendor CORS docs);
  could not capture the header here (egress blocked). A live devtools/curl
  check from an unrestricted network is the recommended final confirmation
  before relying on it either way.

### CONFIRMED — biowasm samtools includes `consensus`

- **Evidence (build recipe):** `biowasm/biowasm` `tools/samtools/compile.sh`
  runs `emmake make samtools` with only `--without-curses` (disables the curses
  TUI `tview`, **not** subcommands) — i.e. a full samtools build, no subcommand
  exclusions. The CDN offers samtools **1.17 / 1.18**, both ≥ 1.16 where
  `consensus` was introduced.
- **Conclusion:** `samtools consensus` is compiled into the WASM module.
- **Confidence:** high for "compiled in." A **runtime smoke test** (load the
  module via Aioli and run `samtools consensus` on a tiny indexed BAM in a
  WebWorker) is still worth doing once, to confirm behaviour end-to-end.

### Correction this verification forces

The earlier lean toward **UCSC** as the reference API (in `genome-api.md`, for
its 0-based coords and `hg19`/`hg38` names) holds **only for a server host**.
For a **static / GitHub Pages** app the reference API must be browser-CORS-safe,
and UCSC apparently is not — so the static-case reference source should be
**Ensembl** (CORS `*`), accepting the `+1` coordinate conversion and the
per-build host split. Summary:

| Host shape | Reference API | Why |
|---|---|---|
| Server (Render etc.) | **UCSC** | server-side call, no CORS needed; nicer coords |
| Static (Pages) | **Ensembl** | CORS `*` confirmed; UCSC needs a proxy |

**Why a server host removes the UCSC CORS problem:** CORS is enforced *only by
the browser* on cross-origin JavaScript requests; server-to-server HTTP is not
subject to it. On a server host the call chain is browser → your backend (same
origin) → UCSC (server-side, no CORS), so UCSC's missing CORS header never comes
into play — exactly how `GnomadApiSource` already calls gnomAD server-side. What
still does *not* work on any host is the **browser** calling UCSC directly;
hosting location doesn't change that, because it's the request *origin* (server
vs. browser), not where the frontend is served, that matters. Caveat for the
server case: all users' UCSC calls share the server's IP, so honour UCSC's
~1 req/s guidance with caching/throttling (it has blocked web apps that didn't),
and Render's free tier cold-starts after idle.

### Still open

- **UCSC CORS — live header capture** from an unrestricted network (devtools or
  `curl -I` with an `Origin` header) to upgrade the UCSC finding from
  "likely/medium-high" to confirmed.
- **Local BAM range-read ergonomics** in Aioli for a user-picked `.bam` + `.bai`
  (File System Access API vs. `File.slice`) on the target browsers.

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
