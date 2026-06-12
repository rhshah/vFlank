# gnomAD GraphQL API — findings, limitations, and design note

Investigation of whether/how to use the public gnomAD API as an alternative
population-frequency masking source (vs downloading per-chromosome VCFs).
Status: **research complete; implementation not started** (awaiting go-ahead).

All facts below were verified against the live endpoint and the gnomAD docs on
2026-06-12.

## Summary / recommendation

Use the **public hosted endpoint** `https://gnomad.broadinstitute.org/api`
(GraphQL over HTTP POST) as an *optional* masking source, selected with
`--pop-source api`. Keep the local-VCF source as the default for reproducible /
bulk / offline (HPC) runs. The API is ideal for **small cohorts with no local
data** (e.g. our 20-variant TP53 test); it is **not** suitable for large cohorts
because of a hard rate limit (below).

No new runtime dependency: use the **standard-library `urllib.request`**.
`requests` is unnecessary; `fastapi` is irrelevant (it builds servers — we are a
client).

## Endpoint and protocol

- URL: `https://gnomad.broadinstitute.org/api`
- Method: HTTP `POST`, `Content-Type: application/json`, body `{"query": "...", "variables": {...}}`
- GraphQL **errors return HTTP 200** with an `errors` array in the body — must be
  checked explicitly, not just by status code.
- Server-side caching present (`x-cached: HIT/MISS` response header).

## Datasets and genome builds (verified)

`DatasetId` enum includes: `gnomad_r4`, `gnomad_r4_non_ukb`, `gnomad_r3` (+subsets),
`gnomad_r2_1`, `gnomad_r2_1_controls`/`_non_neuro`/`_non_cancer`/`_non_topmed`, `exac`.

Mapping for vflank:

| vflank `--genome-build` | `reference_genome` | `dataset` |
|---|---|---|
| `hg19` (GRCh37) | `GRCh37` | `gnomad_r2_1` |
| `hg38` (GRCh38) | `GRCh38` | `gnomad_r4` |

`reference_genome` is the `ReferenceGenomeId` enum (`GRCh37` / `GRCh38`). The
`chrom` argument takes the **bare** name (`"17"`) for both builds — matches our
canonical internal form.

## Query shape (region) and coordinates

```graphql
{
  region(chrom: "17", start: 7579400, stop: 7579540, reference_genome: GRCh37) {
    variants(dataset: gnomad_r2_1) {
      pos ref alt
      exome  { af populations { id ac an } }
      genome { af populations { id ac an } }
    }
  }
}
```

- `region(start, stop)` is **1-based, inclusive**; returned `pos` is 1-based.
  Our masking interface uses 0-based half-open `[s0, e0)` → query `start = s0 + 1`,
  `stop = e0`; returned `pos` values are used directly (they're what `mask_sequence`
  expects).

## Schema facts that matter (introspected)

- `Variant`: `pos, ref, alt, exome, genome, joint, rsids, flags, …`.
- `exome`/`genome` are `VariantSequencingTypeData`: has `af` (Float), `ac`, `an`,
  `populations`, `faf95`, `fafmax`. **`af` IS available at this level.**
- `VariantPopulation`: `id, ac, an, homozygote_count, …` — **no `af` field**; the
  per-population AF must be computed as `ac / an` (guard `an > 0`).
- `joint` (`VariantJointSequencingTypeData`) exists **only in v4** (combined
  exome+genome); it has `ac/an/populations/fafmax` but **no direct `af`** (compute
  `ac/an`). For v2.1.1 `joint` is null.
- v4 uses **`fafmax`** for the filtering-AF max; there is **no `grpmax`** field
  (an early guess that the API corrected).

### AF rule for masking (matches our VCF behaviour of max(AF, AF_grpmax))

For each variant, SNPs only (`len(ref) == 1 and len(alt) == 1`):
```
af_overall = max of (genome.af, exome.af) ignoring None
af_grpmax  = max over (genome.populations + exome.populations) of ac/an, an>0
af = max(af_overall, af_grpmax)
mask pos if af >= af_threshold
```
Verified on the real common SNP TP53 P72R (rs1042522, chr17:7579472 GRCh37):
`exome.af = 0.668`, population-max `= 0.738` — correctly flagged as common.

## Limitations (the important part)

1. **Rate limit: 10 requests / IP / 60 s** (gnomAD policy; README confirms IP-level
   rate limiting in `rate-limiting.ts`). Exceeding it returns errors / blocks the IP.
   This is the dominant constraint.
2. **Per-region variant cap.** A region with too many variants returns
   `"This region has too many variants to display. Select a smaller region."`
   Irrelevant for our ±200 bp windows; relevant if anyone tried whole-gene queries.
3. **Not for bulk.** gnomAD explicitly directs large/programmatic workloads to the
   VCF/Hail-table downloads or the gnomAD Toolbox, not the API.
4. **Reproducibility.** Results are not pinned — gnomAD updates the served data;
   network is required (no offline/HPC compute nodes). For a clinical ddPCR context,
   record the `dataset` id in the run report for provenance.
5. **Schema is "subject to change"** (README) — pin the query and add a schema-shape
   check / clear error if fields go missing.
6. **Coverage gaps.** In coding regions, `genome.af` is often `None` (exome-only
   coverage); the `max-ignoring-None` rule handles this, but it means genome-only
   masking would miss exonic SNPs — use both exome and genome.

### Rate-limit math (why source choice matters)

10 req/min, and our interface calls `get_positions` twice per variant (left + right
flank). Mitigations bring this down:
- **One query per variant window:** pad each region query to cover *both* flanks of
  a variant (`[start-flank-1, end+flank]`) and cache as intervals, so the sibling
  flank call is a cache hit → ~1 request/variant.
- **Region cache:** identical/overlapping variants (e.g. our 20 identical TP53
  rows) collapse to a single request.
- Even so: ~10 *distinct* variant windows/min. ≤~50 distinct variants is comfortable;
  hundreds means minutes of throttled waiting → prefer local VCF.

## Proposed design (implementation pending)

- New `core/popfreq_api.py::GnomadApiSource` implementing the **same duck-typed
  interface** as `GnomadStore`: `get_positions(bare, start0, end0, af_threshold) -> list[int]`.
  Drops in behind `ReferenceFlankSource` with no core changes.
- CLI: `--pop-source {vcf,api}` (default `vcf`). `api` ignores `--pop-vcf-dir`.
- Build→(reference_genome, dataset) mapping as above; `--gnomad-dataset` override
  for subsets (e.g. `gnomad_r2_1_non_cancer`).
- HTTP via `urllib.request`; 30 s timeout; **client-side throttle ≤10/60 s**
  (token bucket); **exponential backoff + retry** on transient/429; raise
  `PopFreqError` with a clear message on hard failure (no silent empty mask).
- Region-interval cache (pad to cover both flanks) to minimise requests.
- Record dataset id + "source=api" in the run report for provenance.
- Up-front guard: if the distinct-variant count is large, warn and suggest
  `--pop-source vcf`.

## Test plan

- **Unit (no network):**
  - `build_to_dataset(genome_build)` mapping.
  - Response→positions parser: given saved JSON fixtures (real API responses for a
    common-SNP region and a rare-only region), assert correct SNP/AF filtering,
    `af = max(overall, ac/an grpmax)`, indels excluded, `None` AF handled.
  - Coordinate conversion (`[s0,e0)` → `start=s0+1, stop=e0`).
  - Throttle/backoff logic with a fake clock and a fake transport (no real HTTP).
  - Error mapping: GraphQL `errors` array and HTTP failure → `PopFreqError`.
- **Fixtures:** commit 2–3 small real JSON responses under `tests/fixtures/gnomad_api/`
  (captured once; <5 KB each) so tests are deterministic and offline.
- **Integration:** the existing `ReferenceFlankSource` masking test, but with a
  stub source object exposing `get_positions` — confirms the API source is
  interchangeable with `GnomadStore` (no real network in CI).
- **Optional live smoke test:** one test hitting the real API, marked
  `@pytest.mark.live` and **skipped by default** (run manually); keeps CI offline
  and polite to the shared resource.

## Docs plan

- `docs/DEVELOPER.md`: add `--pop-source api` to the CLI section and a "new
  frequency source" note (it already documents the duck-typed interface).
- `docs/ARCHITECTURE.md`: add the API source to the masking modes / popfreq seam.
- `ddpcr-conventions` skill: add the AF rule (max of overall + per-pop ac/an) and
  the API-vs-VCF tradeoff so future work reasons about it correctly.
- `README.md`: one line under masking — "no-download option via the gnomAD API for
  small cohorts."
- User-facing guidance everywhere: **API = small/no-download; VCF = bulk/reproducible/HPC.**

## Open decisions for sign-off

1. Default `--pop-source`: keep `vcf` as default (recommended), `api` opt-in.
2. For v4, use `max(genome, exome)` (consistent with v2) or prefer `joint`? Proposed:
   start with `max(genome, exome)` for one code path; revisit `joint` later.
3. Throttle: hard cap at 10/60 s with backoff, or also expose `--api-rate` to let
   power users with a whitelisted IP go faster? Proposed: fixed safe default now.
