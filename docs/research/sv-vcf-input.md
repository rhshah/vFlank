# SV + VCF input — design note (M4 / M4.6)

Captures the breakpoint/strand conventions, the **corrected** fusion-junction
model, and what VCF support requires for both small variants and SVs. Status:
**design; not implemented.** Some items pending user confirmation (flagged ⮕).

## Input formats to support

| Format | Variants | Notes |
|---|---|---|
| MAF | small | already supported (`io/maf.py`) |
| iAnnotateSV / iCallSV TSV | SV | `chr1 pos1 str1 chr2 pos2 str2` (+ optional name/sample) |
| VCF | small **and** SV | one VCF can carry both → classify per record and route |

## Breakpoint strand convention (authoritative: iCallSV)

Per [`iCallSV/dellyVcf2Tab.py`](https://github.com/rhshah/iCallSV/blob/master/iCallSV/dellyVcf2Tab.py):
`str ∈ {0 = top/plus/reference, 1 = bottom/minus/complement}`.

Delly encodes orientation in `INFO.CT = "XtoY"`; map each half `3→0`, `5→1`:

| CT | str1 | str2 | typical |
|----|----|----|---|
| 3to3 | 0 | 0 | inversion |
| 3to5 | 0 | 1 | deletion |
| 5to3 | 1 | 0 | duplication |
| 5to5 | 1 | 1 | inversion |

`chr1=CHROM`, `pos1=POS`, `chr2=INFO.CHR2`, `pos2=INFO.END`. The converter does
**not branch on SVTYPE** — DEL/DUP/INV/INS **and TRA** all go through this single
CT split. TRA (Delly's pre-BND inter-chromosomal type) works because its two
breakpoints take independent orientations, which the two CT halves already
encode. ⮕ *Confirm: Delly always writes TRA's CT as one `XtoY` field (so the
split applies), not two separate fields / paired records.*

Manta/GRIDSS instead use **BND bracket notation** in ALT (`t[p[`, `t]p]`,
`]p]t`, `[p[t`) — a separate parser that must yield the same `(str1, str2)`.
Classify by presence: `INFO.CT` → Delly path; bracketed ALT → BND path.

## Fusion junction model (CORRECTED)

Fused product reads 5'→3' as `partner1 + partner2` (no separator). L =
`fusion_length / 2`. **partner 1 ends at the junction; partner 2 starts at it.**
`pos` is the last base of partner 1 / first base of partner 2 (junction sits
between the two breakpoint bases — confirmed).

| | partner 1 (ends at junction) | partner 2 (starts at junction) |
|---|---|---|
| **str = 0** (+) | `ref[pos−L+1 .. pos]` (+) | `revcomp(ref[pos−L+1 .. pos])` |
| **str = 1** (−) | `revcomp(ref[pos .. pos+L−1])` | `ref[pos .. pos+L−1]` (+) |

Validated on the unambiguous **deletion** case (`CT=3to5` → `str1=0, str2=1`),
whose truth is `plus[..pos1] + plus[pos2..]` with no reverse-complement:
partner 2 must be `+`, not revcomp. An earlier draft reverse-complemented
partner 2 for `str2=1` and was wrong — this table is the fix.

⮕ *Pin with a golden test: one real fusion's `chr1 pos1 str1 chr2 pos2 str2`
plus its expected junction.* (The legacy `.cfg` `partner_dir +/-` is a different,
unverified encoding — do not use it to derive `str`.)

## Small-variant VCF → `Variant`

- Read with `pysam.VariantFile` (bgzip/tabix/headers handled).
- **REF/ALT → start/end** (VCF left-aligns with an anchor base):
  - SNP (`A→T`): `start=end=POS`.
  - insertion (`A→ATG`): anchored at POS (match MAF: `Start=POS, End=POS+1`).
  - deletion (`ATG→A`): deleted span `start=POS+1, end=POS+len(REF)−1`.
  - dedicated tested helper.
- multi-allelic `ALT=A,T` → one `Variant` per ALT.
- header: **always** `Position_REF_ALT`; **plus** gene/HGVS from VEP `CSQ` /
  SnpEff `ANN` when present (blank otherwise).
- **sites-only** — ignore sample IDs / genotypes.

## SV VCF → `Breakpoint` pair → junction

- Delly path: `INFO.CT` (split) + `CHR2` + `END` → `(str1, str2)` as above.
- BND path: parse ALT brackets → same `(str1, str2)`; dedup BND mate via `MATEID`.
- Symbolic `<DEL>/<DUP>/<INV>` with `END`/`SVLEN`: phase 2.

## Variant dedup (small variants)

Today vflank emits **one record per MAF row** → the same variant across N samples
yields N identical records (e.g. the TP53 cohort gave 20). With reference +
population masking the flank is **sample-independent**, so:

⮕ *Confirm:* dedup by `(chrom, start, end, ref, alt)` → **one record per unique
variant**; header = `Position_REF_ALT` + gene/HGVS, **no sample tag** (optionally
`samples=N`). NB: BAM-consensus modes (C/D), being sample-specific, will instead
dedup per `(variant, sample)`.

## Column handling — name-driven, never positional

The simple SV TSV reader picks columns **by header name, not by position**. An
`SvColumns` dataclass holds the logical→header mapping with defaults
(`chr1 pos1 str1 chr2 pos2 str2`, plus optional `name`/`sample`); each is
overridable (mirrors `io/maf.MafColumns` + `--chrom-col` etc.). A row is usable
as long as the header contains the (configured) names — column order is
irrelevant and extra columns are ignored. A missing required column is a clear
error, not a silent positional guess.

## BND parser strategy (phase 2, deferred)

- **Read** VCF with `pysam.VariantFile` (already a dependency; no new dep).
- **Delly path:** read `INFO.CT`/`CHR2`/`END` directly (mirrors iCallSV); no
  library needed.
- **BND-bracket path (Manta/GRIDSS):** implement the ~20-line VCF-spec parser
  ourselves, mapping the 4 forms straight to our `str`, so the mapping stays
  unified with the CT convention and fully unit-tested. Use **vcfpy** /
  **PyVCF** `BreakEnd` and **svtools `vcftobedpe`** as *test oracles*
  (validate against, don't depend on). Crux test: a Manta BND equivalent to a
  deletion must yield `str1=0, str2=1` (= Delly `CT=3to5`).

## Build order

**Phase 1 (now): simple format only.** SV input = the iCallSV/iAnnotateSV TSV,
parsed by column name (above).
1. `core/fusion.py` — `Breakpoint`/`Fusion` + corrected junction model.
2. `io/breakpoints.py` — name-driven TSV reader.
3. `vflank fusion` sub-app + golden test.
4. Small-variant dedup + `Position_REF_ALT` header change.

**Phase 2 (later): VCF.** Small-variant VCF, then Delly-CT / BND SV VCF
(per the BND strategy above).
