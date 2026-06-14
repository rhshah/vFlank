# Plan — true indel-aware consensus (option 3)

Status: **design / not yet implemented.** This is the deferred enhancement
referenced from [bam-consensus.md](bam-consensus.md). Today (option 2) the
consensus is kept **reference-length**: deletions become `N`, and insertion
*anchor* positions are flagged `N` (`insertion_sites()`), but the patient's
actual inserted/deleted bases are not represented. Option 3 makes the consensus
reflect the patient's real template, indels and all.

## Why this matters (and when it doesn't)

A ddPCR primer/probe binds the *physical* patient amplicon, not the reference.
When a patient carries an indel inside a flank, the reference-frame consensus
gives a sequence that does not exist in that sample:

- a 4 bp insertion shifts every downstream base — a probe designed on the
  reference-frame flank is mis-registered past the insertion;
- a deletion masked to `N` throws away usable, designable patient sequence.

Counter-pressure: most robust assays already *avoid* indel-rich regions, and the
reference-frame flank with flagged `N`s is a safe, conservative default. So
option 3 is opt-in, not the new default — it buys accuracy for the assays that
need to design *through* or *next to* a patient indel.

## Hard requirement that makes this non-trivial

Everything downstream assumes a **1:1 reference↔consensus alignment** of length
= window size `N`:

- `apply_lowcov_overlay` indexes `consensus[i]`, `depth[i]`, `reference[i]` by
  the same `i` (genomic position `start + i`);
- gnomAD masking is keyed on **reference** coordinates;
- flank concatenation and `junction_index` (fusions) count reference bases;
- the variant interval `[start, end]` (excluded from flanks) is in reference
  coordinates.

An indel breaks the 1:1 map. The core of option 3 is therefore **representing
the consensus as an explicit reference↔patient alignment**, not a flat string.

## Data model

Replace the flat consensus string with a column list that *is* the alignment:

```python
@dataclass(slots=True)
class ConsensusColumn:
    ref_pos: int | None   # 0-based reference position; None for an inserted column
    base: str             # patient base(s); "" for a pure deletion
    depth: int            # usable depth at the anchoring reference position
    kind: Literal["match", "ins", "del"]

@dataclass(slots=True)
class ConsensusAlignment:
    columns: list[ConsensusColumn]
    start_0based: int     # window start (reference)

    def patient_seq(self) -> str:
        return "".join(c.base for c in self.columns)

    def ref_anchored(self) -> list[ConsensusColumn]:
        return [c for c in self.columns if c.ref_pos is not None]  # len == window size
```

Invariant tests fall straight out of the model:
`len(ref_anchored()) == N`, and
`len(patient_seq()) == N + inserted_bases − deleted_bases`.

## Two candidate engines

**A. Re-align samtools' indel-aware consensus to the reference window.**
Run `samtools consensus` with its default `--show-ins yes --show-del yes`, then
recover the alignment by globally aligning the result back to the reference
window (e.g. `edlib`, which returns a CIGAR cheaply).
*Pros:* reuses the validated htslib caller verbatim.
*Cons:* a new dependency (`edlib`); re-alignment is ambiguous in tandem repeats
(the indel can be placed at several equivalent positions), so the recovered
`ref_pos` map may disagree with the read evidence.

**B. Build the alignment directly from a CIGAR walk (recommended).**
Extend the existing `insertion_sites()` walk (which already proved more reliable
than pysam's pileup stepper on sliced BAMs) into a full per-column consensus:

1. For each reference position in `[start, end)`: tally base calls from spanning
   reads (same primary/MQ/BQ filter as `window_depth`), emit a `match`
   column using `call_fract`/`het_char`; if ≥ `1 − call_fract` of reads carry a
   deletion here, emit a `del` column (`base=""`).
2. At each insertion anchor (≥ threshold, as today), vote the inserted sequence:
   group spanning reads' inserted bases by **length** first, then per-position
   base, and emit one `ins` column carrying the winning inserted string.
3. Apply the low-coverage overlay and gnomAD mask on the `match`/`del` columns
   only (they have a `ref_pos`); `ins` columns are patient-specific and have no
   reference/gnomAD coordinate, so they pass through unmasked.

*Pros:* one engine, explicit columns, no realignment ambiguity, no new
dependency. *Cons:* we reimplement insertion-length voting and must validate it
matches `samtools consensus` on real data.

Recommendation: **B**, validated against samtools' indel-aware output on the
`bam_slice` cohort as an acceptance gate.

## Coordinate-semantics decision

Keep `flank` defined as **N reference bases** (variant anchoring stays in
reference space, gnomAD stays keyed on reference). The *patient* flank that
comes out may be longer or shorter than `N`; both are reported:

- `Flank` (reference span, = `N`) — unchanged meaning;
- `PatientLen` (= `len(patient_seq())`);
- `NInsertedBases`, `NDeletedBases`.

For fusions, `junction_index` becomes a **patient-space** index (`len(patient
partner-1 sequence)`); reverse-complement already operates on the patient string,
so no extra handling — but the genomic-space-before-revcomp masking step must run
on the column model before the string is materialised.

## Surface / flag

Gate behind `--bam-indel {flag,apply}`:

- `flag` (**default**) → today's option-2 behaviour (reference-length, indels
  flagged `N`). No behavior change for existing users.
- `apply` → option 3 (length-changing patient consensus).

## Testing plan

Synthetic BAMs (extend `tests/integration/test_consensus_bam.py`):

- pure insertion → `patient_seq` longer by the insertion; ref-anchored count `N`;
- pure deletion → `patient_seq` shorter; deleted columns dropped;
- het indel (50/50) → resolved per `call_fract`/`het_char`;
- indel adjacent to a common gnomAD SNP → SNP still masked at its `ref_pos`;
- indel straddling the flank/variant boundary → no off-by-one;
- tandem-repeat insertion → placement matches read evidence (engine B's reason
  for being).

Property assertions: `len(ref_anchored) == N`,
`len(patient_seq) == N + inserted − deleted`, and parity with
`samtools consensus --show-ins yes` on cohort windows.

## Phasing

1. **Model + engine B** (`ConsensusAlignment`, CIGAR-walk consensus) with unit
   tests; no wiring yet.
2. **Wire into `ConsensusFlankSource`** behind `--bam-indel apply`; emit the new
   report columns.
3. **Fusions** — patient-space `junction_index`, masking on the column model.
4. **Docs + report** — user guide, update this note to "implemented".

## Risks

- Indel placement ambiguity in repeats (mitigated by engine B using read
  evidence directly).
- Low-coverage overlay at indel boundaries (an `N`-fallback ref column next to an
  `ins` column) — define precedence explicitly: depth gate first, then indel.
- Performance: per-read CIGAR walk over every window — acceptable for the 200 bp
  windows here; revisit if windows grow.
