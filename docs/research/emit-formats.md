# Research & plan — designer-native emit formats (Olivar / Primer3)

Status: **design / not yet implemented.** This is milestone **M4.5-emit** — the
half of the consensus milestone that did not ship in 0.2.0. It is the feature
that makes vflank a *pipeline front-end* rather than a FASTA generator: turning a
masked target into the exact input two downstream designers expect.

## Why this is the priority

Per the project scope (CLAUDE.md): vflank is **not** a primer designer; it must
"add *emit formats* that feed those tools." Today the only output is raw +
`Masked__` FASTA. A user still has to hand-translate that into Olivar's and
Primer3's native inputs — coordinate files, excluded regions, settings blocks —
which is fiddly and error-prone, and is exactly the masking-coordinate
bookkeeping vflank already did internally and then threw away by collapsing to
`N`. Closing this gap is what lets the pipeline actually hand off.

## The core abstraction shift

The single most important design point: **emit needs masked *coordinates*, not
`N`-strings.** Both downstream tools are coordinate-driven (Olivar takes a
variant CSV; Primer3 takes excluded-region offsets). vflank currently computes
the masked positions (`get_positions` / consensus) and then discards them into a
flat `N`-substituted string. The plan introduces a structured intermediate that
all writers consume:

```python
@dataclass(slots=True)
class EmitRecord:
    id: str                       # stable record id (the variant/junction key)
    sequence: str                 # raw (unmasked) sequence
    anchor: tuple[int, int]       # 0-based [start, len) of the variant / junction
                                  #   within `sequence` (the region to design around)
    masked: list[tuple[int, int]] # 0-based (start, len) runs a primer/probe must avoid
    source: str                   # "reference" | "consensus:<sample>"
    meta: dict[str, str]          # gene, HGVSp, rsID(s), sample, AF, …
```

`masked` is the union of common-SNP positions, patient het/low-cov/insertion
sites — coalesced into runs. The current `N`-string output is then just *one*
renderer of this record (`mask_sequence(sequence, masked)`); the FASTA, Olivar,
and Primer3 writers are siblings. This keeps `core/` pure (it produces
`EmitRecord`s) and confines tool-specific knowledge to `io/emit_*.py`.

## Target 1 — Olivar (small-variant amplicons)

Olivar (`treangenlab/Olivar`) designs tiled multiplex amplicon panels and is
variant-aware: it penalises primers that sit over known variable sites. Its
build step consumes a **reference FASTA + a region BED + a variant table**, not
per-target masked sequences. So the Olivar emitter produces:

1. **Target BED** — one interval per variant's flank window
   (`chrom  start  end  name`), giving Olivar the regions to tile.
2. **Variant CSV** — every masked position as a row Olivar can avoid
   (genomic `chrom,pos` + AF/source), reconstructed from `EmitRecord.masked`
   mapped back to genomic coordinates.
3. The reference FASTA is passed through (the user already supplies it via `-r`).

> **Open item:** confirm Olivar's exact variant-CSV column schema and BED
> conventions against the installed Olivar version before coding — Olivar's
> input spec has changed across releases. The emitter must target a pinned
> Olivar version and assert the schema. This is the one place where we cannot
> finalise field names from memory.

## Target 2 — Primer3 (fusion-junction probes, and optionally small variants)

Primer3 reads **Boulder-IO**: `TAG=value` lines, records separated by `=`. The
mapping from an `EmitRecord` is direct and faithful:

```
SEQUENCE_ID={id}
SEQUENCE_TEMPLATE={raw sequence}
SEQUENCE_TARGET={anchor.start},{anchor.len}          # design across the variant/junction
SEQUENCE_EXCLUDED_REGION={s1},{l1} {s2},{l2} …       # the masked runs — primers avoid them
PRIMER_PICK_INTERNAL_OLIGO=1                          # the ddPCR probe (junctions)
=
```

Plus a global settings header (`PRIMER_TASK`, product-size range, Tm/size
constraints) emitted once. Using `SEQUENCE_EXCLUDED_REGION` from the masked
coordinates is **better than feeding `N`s** in the template: Primer3 treats `N`
as a degenerate base (it may still place a primer there), whereas an excluded
region is a hard constraint — which is the masking semantics we actually want.
For fusions, `SEQUENCE_TARGET` is the junction index so the probe spans it.

## Masked vs. unmasked output (the user's question)

Make it explicit and selectable. Add **`--records {both,masked,raw}`** (default
`both`, preserving today's behaviour):

- `both` — raw + masked records (current behaviour; useful for diffing).
- `masked` — only the design substrate (halves output; the common case once you
  trust the masking).
- `raw` — reference-only dump (rare, but a clean reference export).

This applies to the **FASTA** renderer. The emit formats are intrinsically
"masked" (Olivar gets the variant CSV; Primer3 gets excluded regions), so they
are unaffected by `--records` — they always carry the masking as coordinates.

**Coupled semantics fix:** with `--bam`, the `Masked__` record is a *corrected
patient consensus*, not merely masked. The label and `EmitRecord.source` should
distinguish `Masked__` (reference + population mask) from `Consensus__`
(patient). This rename is part of this change, with the old prefix kept as an
accepted alias for one release.

## CLI surface

```
--emit fasta,olivar,primer3      # comma list; default fasta (back-compatible)
--records {both,masked,raw}      # FASTA record selection (default both)
--out-dir DIR                    # emitters write {out-dir}/{format}.{ext}
```

`vflank small` emits FASTA + Olivar + Primer3; `vflank fusion` emits FASTA +
Primer3 (Olivar is amplicon-tiling, not junction probes). Each emitter is a
thin `io/` writer over the shared `EmitRecord` stream.

## Provenance

Every emitted file gets a header comment with the **vflank version**, the
parameter set, and the genome build — the same provenance now printed in the run
report (see [the version-in-report work](#); M-version). A designer artifact you
find six months later must say which vflank built it.

## Module placement

```
core/   produces EmitRecord (pure; no tool knowledge)
io/emit_fasta.py     raw/masked FASTA  (refactor of today's fasta.py renderer)
io/emit_olivar.py    target BED + variant CSV
io/emit_primer3.py   Boulder-IO records + settings header
cli/    --emit / --records wiring; one EmitRecord stream fanned to writers
```

## Testing

- Unit: `EmitRecord` → each writer is a pure string transform — golden-file
  tests (a known variant with two masked SNPs → exact Olivar CSV / Primer3
  block / FASTA).
- Coordinate round-trip: a masked position at genomic `P` must land at the right
  BED row, the right `SEQUENCE_EXCLUDED_REGION` offset, and the right FASTA `N`.
- Integration: emit for the sample MAF, then actually run Olivar/Primer3 on the
  output in CI (behind a marker) to prove the formats parse.
- `--records` matrix: both/masked/raw produce the expected record counts.

## Phasing

1. **`EmitRecord` + refactor FASTA** onto it (no behaviour change; `--records`
   lands here).
2. **Primer3 emitter** (the simpler, fully-specified format) for small + fusion.
3. **Olivar emitter** after confirming its input schema against a pinned version.
4. **Provenance headers + docs + the Olivar/Primer3 round-trip CI marker.**

## Risks / open questions

- **Olivar input schema drift** — pin a version, assert columns (the one
  memory-unsafe spot).
- **Coordinate bases** — BED is 0-based half-open, Primer3 offsets are 0-based
  within the template, MAF is 1-based closed. The `EmitRecord` stores 0-based
  template-relative offsets; genomic round-trips go through the same
  conversions already centralised for flanks.
- **Indels as excluded regions** — a masked insertion anchor is width-0 in the
  reference frame; emit it as a length-1 excluded run around the anchor so
  Primer3 doesn't straddle it.
- **`--records raw` + emit** — raw-only FASTA with an Olivar/Primer3 emit is a
  contradiction (emit needs masking); warn, don't fail.
