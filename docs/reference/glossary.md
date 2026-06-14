# Glossary

Domain and tooling terms used throughout vflank. Acronyms also get hover
tooltips wherever they appear in the docs (via the shared abbreviations list).

## Assay & biology

ddPCR
:   Digital droplet PCR. The partitioned PCR assay that vflank's output is
    designed to feed: a primer/probe must bind a region free of patient/population
    variation, which is why flanks are masked.

Flank
:   The sequence immediately **before** (left) or **after** (right) a variant.
    vflank emits `±N` reference bases around each variant, excluding the variant
    interval itself, as the design substrate for primers/probes.

Masking
:   Replacing a flank base with `N` so downstream design tools avoid it. vflank
    masks common population SNPs (gnomAD) and — with a BAM — patient-specific
    bases that disagree with the reference.

Fusion / junction
:   A chimeric sequence formed by joining two breakpoints from a structural
    variant. vflank builds the reverse-complement-aware junction so a probe can
    span it; the **junction index** marks where partner 2 begins.

Breakpoint
:   One side of a structural variant, given as `chrom:pos:strand`. Strand `0` is
    plus/reference, `1` is minus/complement (the iCallSV convention).

## Variants & coordinates

MAF
:   Mutation Annotation Format. The primary small-variant input. Coordinates are
    **1-based, fully closed** `[start, end]` — unlike pysam's 0-based half-open
    windows, the single most common source of off-by-one bugs (see
    [Architecture](../ARCHITECTURE.md)).

VCF
:   Variant Call Format. Alternative variant input (1-based). Used both for
    small variants and, with Delly `CT`/`BND` records, for structural variants.

SNP
:   A single-nucleotide polymorphism common enough in the population to risk a
    primer/probe mismatch. Only **single-base** population SNPs are masked.

indel
:   An insertion or deletion. In the reference-frame consensus, deletions become
    `N` and insertion anchor sites are flagged `N`; true length-changing indel
    handling is the [option-3 plan](../research/indel-aware-consensus.md).

REF / ALT
:   The reference and alternate alleles. vflank shows the variant literally as
    `[REF/ALT]` between the two flanks in its FASTA output.

## Population frequency (gnomAD)

gnomAD
:   The Genome Aggregation Database — vflank's source for population allele
    frequencies, queried either from a **local VCF** or the **GraphQL API**.

AF / AC / AN
:   Allele Frequency = Allele Count / Allele Number. A flank position is masked
    when its common-SNP AF clears the configured threshold (default `0.001`).

Genome vs exome
:   gnomAD's two callsets. `--pop-data {genome,exome,both}` selects which to
    consult, symmetrically across the VCF and API backends.

Build (hg19/GRCh37, hg38/GRCh38)
:   The genome assembly. AF semantics differ across builds, so the build is
    explicit (`-g`, default `hg19`) and a mismatch warns rather than failing
    silently.

## BAM consensus

BAM
:   A binary alignment file of a sample's reads. With `--bam`/`--bam-map`, vflank
    builds a per-sample patient consensus flank (masking modes C/D).

Consensus
:   The per-base call across a sample's reads, delegated to `samtools consensus`
    (bundled with pysam). Kept reference-length; a pure low-coverage overlay
    falls back to "patient where covered, population where blind."

CIGAR
:   A read's alignment string (M/I/D/N/S…). vflank walks CIGARs directly to flag
    insertions, because pysam's default pileup stepper silently drops reads on
    some sliced BAMs.

MQ / BQ
:   Mapping and base quality thresholds (`--bam-min-mapq`, `--bam-min-baseq`,
    default `20`) gating which reads and bases count toward the consensus.

het / hom
:   Heterozygous / homozygous. A homozygous-ALT call corrects the flank base; a
    heterozygous position becomes `N` (or an IUPAC code with `--bam-het-char
    iupac`).

Low coverage
:   A position below `--bam-min-depth` (default `20`). The default fallback is
    REF + gnomAD masking; `--require-coverage` instead flags the variant.

## Downstream (out of scope for vflank)

Olivar
:   A small-variant amplicon designer. vflank produces its input; it does **not**
    design primers itself.

Primer3
:   A primer/probe designer used for fusion-junction probes. Same boundary:
    vflank emits the junction, Primer3 designs against it.
