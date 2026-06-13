# BAM consensus masking — design note (modes C / D)

Patient-specific flank generation from a sample BAM, for both small variants and
fusions. Status: **design; not implemented.** Decisions taken: full consensus
(not mask-only), `--bam` + `--bam-map`, per-(variant, sample) output.

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
| significant indel evidence | `N` (v1 — see below) |
| **depth < `min_depth`** | see gnomAD layering |

Rationale: a primer over a het site → allele dropout; over a hom-ALT site the
reference would mismatch *every* template, so we use the patient base; low
coverage → we can't confirm, so don't trust it.

## gnomAD layering (mode D)

The principle: **patient-specific where we have coverage; population prior where
we don't.**

- `depth ≥ min_depth` → the BAM is authoritative (table above). gnomAD is ignored
  at that position (if the patient is confidently hom-ALT, a primer matching the
  patient base works even though gnomAD flags the site as common).
- `depth < min_depth` → fall back to gnomAD: `N` if it's a common SNP there, else
  the reference base.

`--bam-lowcov {n,reference,gnomad}` controls the low-coverage base: `n` (safe),
`reference` (complete but unconfirmed), or `gnomad` (the layered default when a
population source is also given). Without a population source, `gnomad` behaves
as `n`.

## Het and indel representation

- `--bam-het-char {N,iupac}` (default `N`). IUPAC (R/Y/S/W/K/M) preserves which
  two bases but most designers treat it as degenerate; `N` is the safe default.
- **Indels (v1 scope):** a flank position with substantial indel evidence is
  masked to `N` — we do **not** reconstruct the indel-shifted sequence (that
  changes flank length / coordinate frame). Indel-aware consensus (à la kindel's
  CIGAR reconciliation) is a deliberate follow-up. The count of indel-masked
  positions is reported (no silent truncation).

## Architecture

The consensus changes the *sequence*, so it is a **FlankSource**, not a mask
source (unlike gnomAD). It reuses the existing `FlankResult` shape:

- `raw` = reference flanks (kept for comparison / the un-`Masked__` record);
- `masked` = the **patient consensus** (hom-ALT corrected, het/low-cov `N`,
  gnomAD-layered at low coverage).

So the 2-records-per-target contract is unchanged — the `Masked__` record simply
becomes the patient-specific assay template.

```
core/consensus.py
  consensus_base(column_stats, reference_base, policy) -> base|'N'   # PURE, unit-tested
  pileup_consensus(bam, chrom, start0, end0, ref_seq, policy)        # pysam pileup loop
BamConsensusSource(bam_path, policy)        # per-sample, lazy-opened + cached
ConsensusFlankSource(reference, bam_source, gnomad=None, flank=...)  # FlankSource (small)
```

Split the pure per-column decision (`consensus_base`) from the pysam pileup
iteration so the policy is testable without a BAM.

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
    seq = pileup_consensus(bam, window, ref_seq, policy)   # patient bases + N
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

- **Synthetic BAM fixtures** (pysam can write SAM/BAM): a het at a known flank
  position → `N`; a hom-ALT → patient base; a clean hom-REF → reference;
  low-coverage → per `--bam-lowcov`; an indel → `N`. The pure `consensus_base`
  gets exhaustive per-policy unit tests with no BAM.
- **Oracle:** validate the position calls against `bcftools mpileup | call` and
  `samtools consensus` on the same synthetic data (run manually / a `@live`-style
  optional test).

## Risks

- **Indels** — the genuinely hard part; v1 masks them (scoped, reported).
- **Build match** — the BAM must be aligned to the same build as `--genome-build`
  / the FASTA; mismatched coordinates give wrong pileups. Add a contig-length
  sanity check against the BAM header.
- **Performance** — per-window pileup is ~ms; cohort cost scales with
  (variants × samples). Cache `(sample, window)`; lazy-open + cache BAM handles.
- **`.bai` required**; clear error if missing.

## Build phases

1. `core/consensus.py` — pure `consensus_base` policy + `pileup_consensus`, with
   synthetic-BAM tests.
2. `BamConsensusSource` + `ConsensusFlankSource`; small-variant integration
   (per-(variant,sample) dedup, sample→BAM mapping, header change).
3. Fusion integration (consensus-before-revcomp in `_segment`).
4. gnomAD layering for low coverage; docs.
