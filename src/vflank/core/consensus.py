"""Patient-specific consensus from a sample BAM (modes C / D).

The base-calling is delegated to ``samtools consensus`` (called in-process via
``pysam.samtools`` — bundled with pysam, no external install). On top of it we
add a **pure** low-coverage overlay: where the BAM is too shallow to trust, fall
back to a gnomAD-informed decision ("patient where covered, population where
blind"). See docs/research/bam-consensus.md.

``apply_lowcov_overlay`` is pure and unit-tested without a BAM; the two pysam
wrappers are exercised against tiny synthetic BAMs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..errors import ConsensusError
from .flanks import FlankResult
from .variant import Variant


@dataclass(slots=True)
class ConsensusPolicy:
    """Thresholds for consensus calling and the low-coverage overlay."""

    min_depth: int = 20          # below this, the BAM call is not trusted
    call_fract: float = 0.9      # fraction of reads to call a (homozygous) base
    het_char: str = "N"          # "N" or "iupac" (ambiguity codes for hets)
    lowcov: str = "gnomad"       # low-coverage base: "n" | "reference" | "gnomad"
    min_baseq: int = 20
    min_mapq: int = 20
    mode: str = "simple"         # samtools consensus model: "simple" | "bayesian"


def _parse_consensus_fasta(text: str) -> str:
    """Join the sequence lines of a single-record FASTA from samtools consensus.

    Deletion placeholders ('*', from ``--show-del yes``) are mapped to 'N' — a
    patient deletion disrupts the flank, so the position is masked.
    """
    seq = "".join(
        line.strip() for line in text.splitlines() if line and not line.startswith(">")
    )
    return seq.upper().replace("*", "N")


def run_samtools_consensus(
    bam_path: str, chrom: str, start_1based: int, stop_1based: int, policy: ConsensusPolicy
) -> str:
    """Return the per-base consensus over a 1-based inclusive region.

    Runs with ``-a`` (emit every position) and ``-d 1`` (call whatever has any
    coverage); the depth-based trust decision is applied separately by
    :func:`apply_lowcov_overlay`, so this returns one base per region position.
    """
    import pysam

    region = f"{chrom}:{start_1based}-{stop_1based}"
    args = [
        "--mode", policy.mode,
        "-a", "-d", "1",
        "-c", str(policy.call_fract),
        "--min-MQ", str(policy.min_mapq),
        "--min-BQ", str(policy.min_baseq),
        # Keep the consensus reference-length (so it aligns to depth/reference and
        # flanks concatenate): drop insertions, mark deletions as '*' -> 'N'.
        # True indel-aware consensus is a deferred enhancement.
        "--show-ins", "no", "--show-del", "yes",
        "-r", region,
    ]
    if policy.het_char == "iupac":
        args.append("--ambig")
    out = pysam.samtools.consensus(*args, bam_path)
    text = out if isinstance(out, str) else "\n".join(out)
    return _parse_consensus_fasta(text)


def window_depth(
    bam_path: str, chrom: str, start_0based: int, end_0based: int, policy: ConsensusPolicy
) -> list[int]:
    """Usable read depth per position over ``[start, end)`` (0-based half-open).

    Counts bases passing base-quality (``min_baseq``) from primary, non-duplicate
    reads with mapping quality ``>= min_mapq`` — an approximation of the depth
    samtools uses, sufficient for the low-coverage gate.
    """
    import pysam

    def keep(read) -> bool:
        return (
            not read.is_unmapped
            and not read.is_secondary
            and not read.is_supplementary
            and not read.is_duplicate
            and read.mapping_quality >= policy.min_mapq
        )

    with pysam.AlignmentFile(bam_path) as af:
        cov = af.count_coverage(
            chrom, start_0based, end_0based,
            quality_threshold=policy.min_baseq, read_callback=keep,
        )
    n = end_0based - start_0based
    return [int(cov[0][i] + cov[1][i] + cov[2][i] + cov[3][i]) for i in range(n)]


def apply_lowcov_overlay(
    consensus_seq: str,
    depth: list[int],
    region_start_0based: int,
    reference_seq: str,
    gnomad_positions: set[int],
    policy: ConsensusPolicy,
) -> str:
    """Override low-coverage positions per the ``lowcov`` policy. Pure.

    All inputs are aligned to the same region: ``consensus_seq[i]``,
    ``depth[i]`` and ``reference_seq[i]`` correspond to 1-based genomic position
    ``region_start_0based + 1 + i``. Where ``depth[i] < min_depth``:
      - ``n``        → ``N``
      - ``reference``→ the reference base
      - ``gnomad``   → ``N`` if that position is a common SNP, else reference.
    Positions with adequate depth keep the samtools call (incl. het ``N``).
    """
    out = list(consensus_seq)
    for i, d in enumerate(depth):
        if d >= policy.min_depth:
            continue
        if policy.lowcov == "n":
            out[i] = "N"
        elif policy.lowcov == "reference":
            out[i] = reference_seq[i].upper()
        else:  # "gnomad"
            pos_1based = region_start_0based + 1 + i
            out[i] = "N" if pos_1based in gnomad_positions else reference_seq[i].upper()
    return "".join(out)


class BamConsensusSource:
    """Per-sample patient consensus from a BAM: samtools engine + low-cov overlay.

    One instance wraps one sample's BAM (lazy-resolved contig notation). Reused
    across that sample's variants; opened handles are cheap.
    """

    def __init__(self, bam_path, policy: ConsensusPolicy | None = None) -> None:
        self.bam_path = str(bam_path)
        self.policy = policy or ConsensusPolicy()
        if not Path(self.bam_path).exists():
            raise ConsensusError(f"BAM not found: {self.bam_path}")
        if not (Path(self.bam_path + ".bai").exists() or Path(self.bam_path + ".csi").exists()):
            raise ConsensusError(
                f"BAM index not found for {self.bam_path} (fix: samtools index {self.bam_path})"
            )
        import pysam

        with pysam.AlignmentFile(self.bam_path) as af:
            self._refs = set(af.references)

    def _contig(self, bare: str) -> str:
        """Resolve a bare chromosome to the BAM's contig name (handles 'chr')."""
        if bare in self._refs:
            return bare
        alt = f"chr{bare}"
        return alt if alt in self._refs else bare

    def consensus(
        self, bare: str, start_0based: int, end_0based: int,
        reference_seq: str, gnomad_positions: set[int],
    ) -> tuple[str, int]:
        """Patient consensus over ``[start, end)``; returns (sequence, n_covered).

        ``n_covered`` is the number of positions at/above ``policy.min_depth``.
        """
        if end_0based <= start_0based:
            return reference_seq.upper(), 0
        contig = self._contig(bare)
        called = run_samtools_consensus(
            self.bam_path, contig, start_0based + 1, end_0based, self.policy
        )
        depth = window_depth(self.bam_path, contig, start_0based, end_0based, self.policy)
        n = end_0based - start_0based
        if len(called) != n or len(depth) != n or len(reference_seq) != n:
            raise ConsensusError(
                f"consensus length mismatch at {bare}:{start_0based}-{end_0based} "
                f"(called={len(called)}, depth={len(depth)}, ref={len(reference_seq)}, want={n})"
            )
        seq = apply_lowcov_overlay(
            called, depth, start_0based, reference_seq, gnomad_positions, self.policy
        )
        return seq, sum(1 for d in depth if d >= self.policy.min_depth)

    def close(self) -> None:  # no persistent handle; symmetry with other sources
        pass


class ConsensusFlankSource:
    """FlankSource producing patient-consensus flanks (raw=reference, masked=consensus).

    gnomAD (if given) provides the low-coverage fallback decision per the policy.
    """

    def __init__(self, reference, bam_source: BamConsensusSource, gnomad=None,
                 *, flank: int = 200, af_threshold: float = 0.001):
        self.reference = reference
        self.bam = bam_source
        self.gnomad = gnomad
        self.flank = flank
        self.af_threshold = af_threshold

    def _gnomad_positions(self, bare: str, start_0based: int, end_0based: int) -> set[int]:
        if self.gnomad is None:
            return set()
        return set(self.gnomad.get_positions(bare, start_0based, end_0based, self.af_threshold))

    def fetch(self, variant: Variant) -> FlankResult:
        left_start_0 = max(0, variant.start - self.flank - 1)
        left_end_0 = variant.start - 1
        right_start_0 = variant.end
        right_end_0 = variant.end + self.flank

        left = self.reference.fetch(variant.chrom, left_start_0, left_end_0)
        right = self.reference.fetch(variant.chrom, right_start_0, right_end_0)
        masked_left, cov_left = self.bam.consensus(
            variant.chrom, left_start_0, left_end_0, left,
            self._gnomad_positions(variant.chrom, left_start_0, left_end_0),
        )
        masked_right, cov_right = self.bam.consensus(
            variant.chrom, right_start_0, right_end_0, right,
            self._gnomad_positions(variant.chrom, right_start_0, right_end_0),
        )
        return FlankResult(
            left.upper(), right.upper(), masked_left, masked_right,
            covered=cov_left + cov_right, total=len(left) + len(right),
        )
