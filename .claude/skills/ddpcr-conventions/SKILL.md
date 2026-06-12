---
name: ddpcr-conventions
description: ddPCR assay-design domain reference for the vflank repo — genomic coordinate systems (1-based MAF vs 0-based pysam), MAF column semantics, gnomAD allele-frequency masking, flank masking modes A–D, heterozygous/consensus handling, and fusion-junction logic. Use when editing or reasoning about flank extraction, coordinate math, SNP/population-frequency masking, BAM consensus, or fusion breakpoints in this codebase.
---

# ddPCR conventions & domain reference

Background knowledge for working on vflank correctly. The code-level conventions
live in `CLAUDE.md`; this is the *why* and the biology that prevents silent,
expensive mistakes.

## What ddPCR needs from us, and why masking matters

Droplet Digital PCR quantifies a target by partitioning DNA into thousands of
droplets and amplifying with a primer pair + a sequence-specific probe. The
assay only works if the **primers and probe anneal cleanly to the template that
is actually present in the patient.**

The failure mode vflank exists to prevent: if a primer or probe sits over a
**polymorphic position** (a common germline SNP, or a patient-private variant),
then in patients carrying the alternate allele the oligo mismatches → reduced
binding, allele dropout, or assay failure. So we **mask** positions that vary,
turning them into `N` (or an IUPAC code) so the downstream designer avoids
placing an oligo 3′ end there. We mask the **flanks** (primer-landing territory);
the variant of interest itself is always shown literally as `[REF/ALT]` — that's
the target, not something to hide.

Default flank is ±200 bp, matching ddPCR amplicon scale (~60–200 bp products,
designed from a wider candidate window).

## Coordinate systems — get this wrong and the sequence is silently wrong

- **MAF**: 1-based, fully-closed `[start, end]`. A SNP has `start == end`.
- **pysam** (`FastaFile.fetch`, tabix `fetch`): 0-based, half-open `[start, end)`.
- **VCF** `POS`: 1-based.

Flank extraction (`core/flanks.ReferenceFlankSource`):
```
left  = fetch(chrom, max(0, start - flank - 1), start - 1)   # bases before variant
right = fetch(chrom, end, end + flank)                        # bases after variant
```
Masking a 1-based VCF position into a flank string:
```
idx = pos - region_start_0based - 1
```
**pysam never raises on an over-run** — it returns a truncated string at contig
ends. A flank shorter than requested is a real condition; report it, never drop
the record silently.

## MAF column semantics (TCGA/MSK)

Required: `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`,
`Tumor_Seq_Allele2` (the somatic alt — *Allele2*, not Allele1). Metadata used in
headers: `Hugo_Symbol`, `HGVSp_Short`, `HGVSc`, `Tumor_Sample_Barcode`.

Allele encoding: `-` denotes an empty allele (insertion has `Reference_Allele = -`;
deletion has `Tumor_Seq_Allele2 = -`). Chromosome may arrive as `7`, `chr7`,
`Chr7`, numeric `23/24/25` (X/Y/MT), or `M/chrM`; `core/chrom.normalise_chrom`
canonicalises all of these to the bare form.

## gnomAD allele-frequency masking

We mask a flank position if a gnomAD variant there is a **single-base
substitution** (REF length 1, every ALT length 1) whose **max AF ≥ threshold**
(default 0.001 = 0.1%). Indels are deliberately *not* masked (conservative —
avoids over-masking primer territory).

AF fields differ by build/release:
- **gnomAD v4.1 (hg38)** exposes both `AF` and `AF_grpmax`; we take the max.
- **gnomAD v2.1.1 (hg19)** has `AF` but **no `AF_grpmax`** — only `AF` is used.
- Guard against `.` and `NaN` AF tokens (the parser rejects `NaN` via `f == f`).

Files are per-chromosome bgzipped VCFs + tabix index; resolution is by known
filename patterns (`core/popfreq.GNOMAD_PATTERNS`). A contig legitimately absent
from a VCF (e.g. MT) is logged at DEBUG, not failed.

## Flank source modes (the `FlankSource` strategy)

| Mode | Inputs | Flank source | Masking |
|------|--------|--------------|---------|
| A | MAF + FASTA | reference | none |
| B | + gnomAD | reference | common SNPs → N |
| C | + sample BAM | patient consensus, reference fallback at low depth | het / low-confidence → N or IUPAC |
| D | all | patient consensus | gnomAD ∪ observed-het |

A and B are implemented (`ReferenceFlankSource`). C/D are the differentiator:
patient consensus catches **private/rare** variants gnomAD never sees — the ones
that silently break a primer for one specific patient.

### Consensus (modes C/D) plan
- Build from BAM via `bcftools mpileup → call → consensus --iupac-codes`, or
  pysam pileup. Validate against `samtools consensus` as an oracle.
- A position is **heterozygous** (→ mask) if the second-most-common allele
  exceeds ~25–30% with sufficient depth; below a min-depth threshold, fall back
  to the reference base and flag it.
- Indels in the pileup shift coordinates — the hard part. `kindel` is the
  reference implementation for CIGAR-described indel reconciliation, but it is
  haploid/clonal (viral) and does **not** flag heterozygosity, so its majority
  call must be augmented with diploid het detection.

## Fusion / structural-variant path (SV)

A gene fusion is defined by two breakpoints on (possibly) different chromosomes.
The chimeric **junction sequence** is the only fusion-specific template, so the
ddPCR probe must **span the junction**. Assembly: fetch each partner's flank,
orient by transcript strand (reverse-complement the `-` partner), concatenate
across the breakpoint — **without inserting any non-ACGT separator into the
sequence** (the legacy script's `"-"` join corrupted the junction base; do not
reproduce it). Probe design over the junction is delegated to Primer3, honouring
the config's GC/Tm/length (the legacy script ignored its own GC settings).

Config (`config_ES_CTDNA_03.cfg`) is hg19/GRCh37 and describes one fusion per
run (e.g. EWSR1–WT1 for Ewing/DSRCT). The breakpoint orientation semantics in
the legacy script are ambiguous and must be recovered from a domain expert before
re-implementation.

## Downstream emit formats (where vflank stops)

- **Olivar** (small-variant amplicons) consumes a FASTA + a SNP CSV with columns
  `START, STOP, FREQ` (1-based). vflank's gnomAD scan and BAM-het detection both
  produce exactly this — `--emit-olivar` is the integration seam. The `FREQ`
  column is where patient-specific risk (from a BAM) gets injected.
- **Primer3** handles the fusion-junction probe.

vflank produces inputs for these tools; it does not design primers itself.
