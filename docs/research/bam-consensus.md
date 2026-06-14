# BAM consensus masking — design note (modes C / D)

Patient-specific flank generation from a sample BAM, for both small variants and
fusions. Status: **implemented** (small + fusion). User guide:
[user-guide/consensus.md](../user-guide/consensus.md). Decisions taken: full
consensus (not mask-only), `--bam` + `--bam-map`, per-(variant, sample) output.

## Goal

gnomAD masking is a *population* prior. With the patient's own reads we can do
better: build the flank from the **patient's actual sequence** so a primer/probe
matches the real template, and `N` only where the patient is genuinely
ambiguous. This catches **private/rare** variants gnomAD never saw.

- **Mode C** = BAM consensus only.
- **Mode D** = BAM consensus, with gnomAD as the fallback where the BAM is blind.

## Per-position consensus policy

For each flank position, pile up reads passing `min_mapq` / `min_baseq`, with
`ignore_overlaps=True` (don't double-count overlapping mates). Then:

| Condition (depth ≥ `min_depth`) | Output base |
|---|---|
| homozygous-REF (ref frac ≥ `consensus_fraction`, e.g. 0.9) | reference base |
| homozygous-ALT (one alt ≥ `consensus_fraction`) | **patient ALT base** (consensus correction) |
| heterozygous (minor allele ≥ `het_fraction`, e.g. 0.25) | `N` (or IUPAC) |
| patient indel | reflected in the consensus (samtools is indel-aware) |
| **depth < `min_depth`** | see gnomAD layering |

This policy is realised through `samtools consensus` options (`--call-fract` →
hom calls, default/`--ambig` → het, `--min-depth`) plus our low-coverage overlay
(below). Rationale: a primer over a het site → allele dropout; over a hom-ALT
site the reference would mismatch *every* template, so we use the patient base;
low coverage → we can't confirm, so don't trust it.

## gnomAD layering (mode D)

The principle: **patient-specific where we have coverage; population prior where
we don't.**

- `depth ≥ min_depth` → the BAM is authoritative (table above). gnomAD is ignored
  at that position (if the patient is confidently hom-ALT, a primer matching the
  patient base works even though gnomAD flags the site as common).
- `depth < min_depth` → fall back to gnomAD: `N` if it's a common SNP there, else
  the reference base.

`--bam-lowcov {n,reference,gnomad}` controls the low-coverage base: `n` (mask),
`reference`, or `gnomad` (**default**: the reference base, with gnomAD N-masking
common SNPs when a population source is given). So below `min_depth` (default
20×) we fall back to **REF + gnomAD** — an uncovered variant just behaves like a
normal no-BAM run (not all-N). Without a population source, `gnomad` is plain
`reference`.

## Het and indel representation

- `--bam-het-char {N,iupac}` (default `N`) → maps to `samtools consensus`
  default (`N`) vs `--ambig` (IUPAC R/Y/S/W/K/M). `N` is the safe default for
  primer design.
- **Indels: reference-frame (v1).** `samtools consensus` is indel-aware, but its
  indel output changes the sequence length, which breaks the position-by-position
  overlay and flank concatenation. So we run it with `--show-ins no --show-del
  yes` to keep the consensus **reference-length**. Both indel kinds are flagged,
  not lost:
  - **Deletions** become `*` → `N` directly in the engine output (a patient
    deletion disrupts the flank, so the position is masked).
  - **Insertions** would otherwise be dropped silently by `--show-ins no`. To
    avoid a silent failure we independently scan read CIGARs
    (`insertion_sites()`, same primary/MQ filter as the depth array) and **mask
    the anchor base — the reference position immediately before the insertion —
    to `N`** when the insertion is supported by ≥ `1 − call_fract` of spanning
    reads at ≥ `min_depth`. The count is surfaced as `NInserted` per variant and
    an `Insertion sites` line in the run summary.

  Reads are scanned via `fetch` + CIGAR walk rather than pileup: pysam's default
  pileup stepper silently drops reads in some BAMs (e.g. subset/sliced BAMs),
  which would under-count insertions. True indel-aware (length-changing)
  consensus — keeping the patient's actual inserted bases — is a deferred
  enhancement (see [option 3 plan](indel-aware-consensus.md)).

## Engine — `samtools consensus` via pysam (hybrid)

The base-calling is delegated to **`samtools consensus`**, called in-process via
`pysam.samtools.consensus(...)` — bundled with pysam (verified; no extra system
dependency; `bcftools` is *not* bundled in this build). It is the validated
htslib implementation and is indel-aware, so we do not reimplement consensus.

We add a thin, **pure** overlay for the two things samtools doesn't do:

1. **gnomAD low-coverage layering** — get a depth array over the window
   (`pysam.AlignmentFile.count_coverage`, in-process), and at positions below
   `--bam-min-depth` replace the call with the gnomAD decision (`N` if a common
   SNP there, else the reference base). "Patient where covered, population where
   blind."
2. **Variant-of-interest exclusion** — already handled (flank windows are
   strictly outside `[start, end]`).

```
core/consensus.py
  run_samtools_consensus(bam, chrom, start, stop, opts) -> str   # impure (pysam.samtools)
  window_depth(bam, chrom, start0, end0) -> list[int]            # impure (count_coverage)
  apply_lowcov_overlay(seq, depth, region_start0, ref_seq,       # PURE, unit-tested
                       gnomad_positions, min_depth) -> str
BamConsensusSource(bam_path, policy, engine="samtools")  # per-sample, lazy-open + cache
ConsensusFlankSource(reference, bam_source, gnomad=None, flank=...)  # FlankSource (small)
```

The consensus changes the *sequence*, so `ConsensusFlankSource` is a
**FlankSource** (not a mask source). It reuses the existing `FlankResult`:
`raw` = reference flanks; `masked` = the patient consensus. The
2-records-per-target contract is unchanged — `Masked__` simply becomes the
patient-specific assay template.

**Pluggable engine** (`--bam-engine {samtools,pileup}`, default `samtools`): the
default delegates to samtools; a future pure-pysam-pileup engine can be added
behind the same `BamConsensusSource` interface for in-process speed — exactly
the pattern used for the gnomAD `vcf`/`api` backends. The overlay is engine-independent.

## Performance & parallelism

`samtools consensus` is a subprocess per region, so cost scales with
(variants × samples). This is **not** optimised inside vflank — the tool stays a
clean per-unit processor and **Nextflow fans the work out** (per-sample,
per-chromosome, or per-region) at the pipeline layer. The swappable `pileup`
engine remains an option if single-process throughput ever matters.

## Sample → BAM mapping

- `--bam SAMPLE.bam` — single sample (require a one-sample MAF, or combine with
  `--samples`).
- `--bam-map samples.tsv` — `sample<TAB>bam_path` for cohorts (explicit/auditable).
- A sample with no BAM → **warn**, fall back to reference + gnomAD (mode B) for
  that sample (never silent).

## Output granularity (the data-model change)

BAM consensus is **sample-specific**, so dedup changes:

- **With any `--bam*`:** dedup key = `(CHR_POS_REF_ALT, sample)` → one record per
  (variant, sample); the **sample returns to the header**
  (`>{SAMPLE}__{GENE}__{HGVSp}__{HGVSc}__{CHR}_{POS}_{REF}_{ALT}`). Applied to all
  records for a consistent granularity, even BAM-less samples (which get mode B).
- **Without `--bam*`:** unchanged — dedup by `CHR_POS_REF_ALT`, no sample tag.

## End-to-end — small variants

```
for each (variant, sample):
    bam = bam_for(sample)                     # --bam or --bam-map; may be None
    raw    = reference flanks
    masked = ConsensusFlankSource(reference, bam, gnomad).fetch(variant)
             # per-position policy + gnomAD layering at low coverage
emit (raw, masked); dedup key (CHR_POS_REF_ALT, sample); sample in header
```

## End-to-end — fusions

Each partner flank is an ordinary reference window (reads spanning the *junction*
won't align — that's fusion evidence — but the partner flanks pile up normally).
The consensus is built in **genomic (plus-strand) space, then reverse-complemented**
— our existing trick: `revcomp` of a consensus segment is just the complement
(real bases complement; `N`→`N`), so consensus-before-revcomp is correct.

```
_segment(reference, bp, flank, *, donor, bam_source, gnomad):
    window = genomic flank window
    seq = bam_source.consensus(window, ref_seq, gnomad)    # samtools + lowcov overlay
    return revcomp(seq) if rc else seq
junction = partner1 + partner2 ; one record per (fusion, sample)
```

The breakpoint TSV's optional `sample` column drives `bam_for(sample)`.

## CLI surface

```
vflank small run MAF -r REF -g hg19 \
    --bam-map samples.tsv \            # or --bam single.bam
    --bam-min-depth 20 --bam-het-fraction 0.25 \
    --bam-consensus-fraction 0.9 --bam-min-baseq 20 --bam-min-mapq 20 \
    --bam-het-char N --bam-lowcov gnomad \
    --pop-source api                   # gnomAD fallback for low-coverage
# fusion: identical --bam* flags; uses the TSV sample column
```

## Testing & validation

- **Synthetic BAM fixtures** (pysam writes SAM/BAM in-process): het → `N`/IUPAC;
  hom-ALT → patient base; hom-REF → reference; low-coverage → per `--bam-lowcov`;
  indel → reflected (not masked). End-to-end through `pysam.samtools.consensus`.
- The **pure `apply_lowcov_overlay`** gets exhaustive unit tests with no BAM.
- **Oracle:** the engine *is* `samtools consensus`, so it is self-validating for
  the calling; the overlay is what we test. `bcftools` is unavailable in-process
  (not bundled), so any bcftools cross-check is a manual/optional step.

## Risks

- **Build match** — the BAM must be aligned to the same build as `--genome-build`
  / the FASTA; mismatched coordinates give wrong calls. Sanity-check contig
  lengths against the BAM header.
- **Performance** — subprocess per region; not optimised in vflank (Nextflow fans
  out; `pileup` engine is the in-process fallback). Lazy-open + cache BAM handles;
  cache `(sample, window)`.
- **`.bai` required**; clear error if missing.
- **samtools availability** — bundled with pysam, so no external install; pin a
  pysam floor that includes `samtools consensus` (samtools ≥ 1.16).

## Build phases

1. `core/consensus.py` — `run_samtools_consensus` + `window_depth` +
   pure `apply_lowcov_overlay`, with synthetic-BAM tests.
2. `BamConsensusSource` + `ConsensusFlankSource`; small-variant integration
   (per-(variant,sample) dedup, `--bam`/`--bam-map`, header change).
3. Fusion integration (consensus-before-revcomp in `_segment`).
4. gnomAD low-coverage overlay polish + docs; `--bam-engine` knob (samtools
   default; `pileup` engine deferred).
