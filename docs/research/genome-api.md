# Reference-sequence API — findings and design note

Investigation of whether/how to fetch the **reference sequence** from a public
API instead of a local FASTA, so a hosted/no-download run needs no multi-GB
reference on disk. Status: **implemented** — `core/reference_api.ReferenceApiSource`,
selected via `--ref-source api`. Symmetric with the gnomAD story
(`--pop-source vcf|api`): a local `ReferenceFasta` (default) and an API source
behind the same duck-typed surface `ReferenceFlankSource` depends on.

## Summary / recommendation

Use the **UCSC getData/sequence endpoint** (`https://api.genome.ucsc.edu/getData/sequence`)
as an *optional* reference source, selected with `--ref-source api`. Keep the
local FASTA as the default for reproducible / bulk / offline (HPC) runs. The API
is ideal for the **hosted single-variant / small** case where shipping a FASTA
is impractical (e.g. a free-tier container).

UCSC was chosen over Ensembl REST for two correctness-ergonomics reasons:

| | **UCSC** (chosen) | Ensembl REST |
|---|---|---|
| Build id | `hg19` / `hg38` — *our `--genome-build` values verbatim* | GRCh38 main host; GRCh37 on `grch37.rest.ensembl.org` |
| Coordinates | **0-based half-open** — identical to pysam / our flank math | 1-based inclusive (needs +1) |
| Contig | `chr`-prefixed (`chr17`, `chrM`) | bare (`17`) |
| Browser CORS | not browser-safe (server-side only) | `Access-Control-Allow-Origin: *` |

Because vflank calls the reference **server-side**, UCSC's lack of browser CORS
is irrelevant; the 0-based-coords match keeps the flank math in `flanks.py`
unchanged (the #1 coordinate-bug risk). Ensembl would only matter for a future
**static/browser** build (no server) — deferred with the rest of that path; see
[web-app-and-hosting.md](web-app-and-hosting.md).

No new runtime dependency: standard-library `urllib.request`, mirroring
`popfreq_api`.

## Endpoint and protocol

- URL: `https://api.genome.ucsc.edu/getData/sequence`
- Method: HTTP `GET`, query params `genome`, `chrom`, `start`, `end`.
- `start`/`end` are **0-based half-open** — our internal `[s0, e0)` is passed
  through unchanged (no conversion, unlike the gnomAD region query).
- Response: JSON `{ "genome", "chrom", "start", "end", "dna", ... }`. An error
  payload carries an `"error"` string.
- Contig naming: UCSC uses `chr`-prefixed names; the mitochondrion is **`chrM`**
  (not `chrMT`). `ucsc_contig()` maps our bare `MT`/`M` → `chrM`, else `chr{bare}`.
  We cannot cheaply probe the contig list over the API, so we apply UCSC's known
  convention rather than auto-detecting (as `ReferenceFasta.contig` does).

## Masking-rule parity / coordinate parity

`ReferenceApiSource.fetch(bare, s0, e0)` returns the bases for `[s0, e0)` — the
exact contract of `ReferenceFasta.fetch`, so `ReferenceFlankSource` is unchanged.
A **shorter-than-requested** sequence (contig end) is returned as-is, *not*
raised — identical to a pysam short read — and the existing CLI flank-truncation
check reports it. The only hard-error paths are an `error` payload, a missing
`dna` field, non-JSON, or exhausted retries → `ReferenceError`.

## Build handling (a deliberate difference from the FASTA source)

`ReferenceFasta` fingerprints the build by chr1 length to catch hg19/hg38
mix-ups. With an API reference there is no local sequence to fingerprint, so
`ReferenceApiSource.check_build` returns `None` and **trusts** the requested
build (passed straight through as the UCSC `genome`). A wrong build still
surfaces — as a UCSC error or wrong-length window — rather than as silent wrong
sequence, but the proactive chr1-length warning is not available. The build
selection is therefore the user's responsibility (in the hosted UI, the
GRCh37/GRCh38 form control).

## Limitations

1. **Rate limit.** UCSC asks for light use (~1 request/second) and has blocked
   web apps that exceeded it. Mitigated by a client-side throttle (`min_interval=1.0`)
   + a per-window cache; fine for single-variant / small, not for bulk.
2. **Not for bulk.** Hundreds of variants → minutes of throttled waiting; use a
   local FASTA.
3. **Network required** (no offline/HPC compute nodes).
4. **Shared-IP throttle when hosted.** All users' calls leave one server IP, so
   the single-instance throttle + cache is also the rate limiter (see
   web-app-and-hosting.md, Render note).

## Design (as implemented)

- `core/reference_api.py::ReferenceApiSource` — same surface the CLIs use on a
  reference: `fetch(bare, s0, e0) -> str`, `check_build`, `has_chr`, `close`,
  plus `request_count` for monitoring. Drops in behind `ReferenceFlankSource`.
- Pure, unit-tested kernels: `build_url`, `parse_sequence_response`,
  `genome_for_build`, `ucsc_contig`. HTTP/throttle/clock injected for offline tests.
- Throttle (≤1 req/s) + linear backoff/retry on transient (network/timeout/429/5xx),
  window cache, typed `ReferenceError` on hard failure — no silent empty result.
- CLI: `cli/_reference.make_reference_source` + `--ref-source {file,api}`
  (default `file`); `--ref-genome` is now optional and required only for `file`
  (clear error otherwise — no silent fallback). Mirrors `_masking`.
- Monitoring: `Reference API req` in the run summary and `reference_api_requests`
  in the `--report` stats.

## Open items (verify before relying on the API in anger)

- **Live response-shape check.** This was implemented from UCSC's *documented*
  response shape; the session that wrote it had no egress to `api.genome.ucsc.edu`
  (proxy-blocked), so the JSON shape and the `chrM` mito name are **not yet
  live-verified**. Run the `@pytest.mark.live` smoke test
  (`VFLANK_LIVE_API=1 pytest tests/unit/test_reference_api.py -k live`) from an
  unrestricted network and confirm: `dna` field present, coordinates 0-based,
  `chrM` for mitochondria, soft-masked (lowercase) bases handled by the existing
  uppercasing.
- **UCSC rate-limit behaviour** under the hosted shared IP (confirm the ~1 req/s
  throttle is sufficient, or add caching headroom).

## Test plan (as implemented)

- **Unit (no network)** — `tests/unit/test_reference_api.py`: `genome_for_build`,
  `ucsc_contig` (incl. `chrM`), `build_url` coords, `parse_sequence_response`
  (ok / error payload / missing `dna` / short sequence), fetch+cache+throttle,
  retry-then-success, exhausted-retry → `ReferenceError`, `check_build` trust.
- **Integration** — `tests/integration/test_small_pipeline.py`: `small run
  --ref-source api` via `CliRunner` with an injected fake transport (no network,
  no FASTA), asserting flank output + the monitoring line; and `--ref-source
  file` without `--ref-genome` → clean error.
- **Live smoke test** — one test hitting the real UCSC API, `@pytest.mark.live`,
  skipped by default.
