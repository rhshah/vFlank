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
    """Join the sequence lines of a single-record FASTA from samtools consensus."""
    return "".join(
        line.strip() for line in text.splitlines() if line and not line.startswith(">")
    ).upper()


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
