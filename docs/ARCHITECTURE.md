# vflank — Architecture & Roadmap

## What vflank is (and is not)

vflank is the **variant-aware, optionally patient-specific, masked-flank
front-end** of a ddPCR assay-design pipeline. It extracts the sequence flanking
each variant and masks positions that would compromise a primer/probe, then
emits clean targets in designer-native formats.

It is **not** a primer designer. Design is delegated downstream:

| Path | Designer | Why |
|------|----------|-----|
| Small variants (SNP/indel) | **Olivar** (SADDLE, risk-array) | Built for variant-aware multiplex/tiled amplicons. |
| Fusion junction | **Primer3** | Single junction-spanning hybridisation probe. |

Olivar is **GPL-3.0** and pulls in BLAST/MAFFT/NumPy<2, so it is invoked
**out-of-process** (CLI / container), keeping vflank permissively licensed and
the dependency graph clean. Orchestration across samples/tools is a **Nextflow**
layer that wraps the stable CLIs in per-process containers — added last, never a
prerequisite for the standalone CLI.

```
                 ┌── small variants ──► Olivar  ─┐
vflank front-end ─┤   (FASTA + SNP CSV)           ├─► ddPCR assays
 (flank + mask)   └── fusion junction ─► Primer3 ─┘
        ▲
   ref FASTA  ·  gnomAD VCFs  ·  (optional) sample BAM
```

## Flank source strategy (modes)

A `FlankSource` decides *where each flank base comes from*:

| Mode | Inputs | Source | Masking |
|------|--------|--------|---------|
| A Reference | MAF + FASTA | reference | none |
| B Reference + pop-mask | + gnomAD | reference | common SNPs → N |
| C Consensus | + BAM | patient consensus, ref fallback | het / low-cov → N/IUPAC |
| D Consensus + pop-mask | all | patient consensus | gnomAD ∪ observed-het |

Mode C/D is the differentiator: patient consensus catches **private/rare**
variants gnomAD never sees — the ones that silently break a primer for one
patient. Consensus will be built on `bcftools mpileup→call→consensus
--iupac-codes` (or pysam pileup), validated against `samtools consensus`, with
kindel's indel reconciliation as a reference for the hard CIGAR cases.

Implemented today: `ReferenceFlankSource` (modes A + B). The population-mask
backend is pluggable behind a duck-typed `get_positions` interface, with two
interchangeable implementations selected by `--pop-source`:
`GnomadStore` (local VCFs) and `GnomadApiSource` (gnomAD GraphQL API, no
download, rate-limited). Both honour `--pop-data {genome,exome,both}` (union for
`both`). See [research/gnomad-api.md](research/gnomad-api.md).

## Module map

```
src/vflank/
├── core/
│   ├── chrom.py     notation detect/normalise (pure)
│   ├── variant.py   Variant dataclass + validation (pure)
│   ├── flanks.py    FlankSource protocol, ReferenceFlankSource, mask_sequence
│   ├── popfreq.py   gnomAD VCF resolve + parse_common_snp_positions (pure) + GnomadStore
│   └── popfreq_api.py  gnomAD GraphQL API source (GnomadApiSource) + pure parser
├── io/
│   ├── maf.py       load/remap/validate, row → Variant
│   ├── reference.py ReferenceFasta + genome-build guard
│   └── fasta.py     header sanitise + record format/write
├── cli/
│   ├── app.py       root Typer, global -v/-q/--debug, version
│   └── small.py     run · inspect · list-vcf
├── logging.py       Rich console + logger
└── errors.py        VflankError hierarchy
```

The hot kernels (`parse_common_snp_positions`, future consensus pileup) are pure
functions over plain iterables — the natural seam to later accelerate with Rust
(rust-htslib for consensus parity with samtools; noodles for pure-Rust gnomAD
scanning). Lock correctness in Python first; port the proven bottleneck only.

## Origin

vflank began as two scripts: `get_flanking_sequence.py` (small variants, fully
ported in M2 with the output format preserved) and `design_fusion_primers.py` +
`config_ES_CTDNA_03.cfg` (a legacy Python-2 fusion script that didn't run and
embedded a junction-corrupting `"-"`). Both have been ported / re-implemented and
removed; the conventions recovered from them (the corrected fusion-junction
model, iCallSV strand mapping) are captured in
[research/sv-vcf-input.md](research/sv-vcf-input.md).

## Roadmap

| M | Deliverable | State |
|---|-------------|-------|
| M1 | Scaffold (src-layout, pyproject/Hatch, Typer root, logging) | ✅ done |
| M2 | Small-variant port (behaviour-preserving, testable) | ✅ done |
| M3 | Tests + provenance (unit + integration, build guard) | 🟡 in progress |
| M4 | Fusion rewrite (Python 3 + Pydantic config) | ⬜ |
| M4.5 | Consensus (modes C/D) + `--emit-olivar` / `--emit-primer3` | ⬜ |
| M5 | MkDocs (Material + mkdocstrings) | ⬜ |
| M6 | Release (PyPI via OIDC; pysam wheel caveats) | ⬜ |
| M7 | Nextflow / nf-core pipeline (containerised) | ⬜ |

### Top risks
1. Fusion script intent is undocumented/buggy — recover before M4.
2. Reference build mismatch → silent wrong sequence — guarded in `reference.py`.
3. pysam wheel fragility — CI matrix Linux+macOS; document WSL/conda.
4. Indel coordinate edge cases — dedicated tests in M3.
5. Olivar GPL/dep isolation — enforced by the out-of-process boundary.
