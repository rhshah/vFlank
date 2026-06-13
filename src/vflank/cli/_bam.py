"""CLI helpers for BAM consensus: policy construction and sample->BAM resolution.

Shared by ``small`` and (later) ``fusion``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from ..core.consensus import ConsensusPolicy
from ..errors import VflankError


def build_consensus_policy(
    min_depth: int, call_fract: float, het_char: str, lowcov: str,
    min_baseq: int, min_mapq: int,
) -> ConsensusPolicy:
    if het_char not in ("N", "iupac"):
        raise VflankError(f"--bam-het-char must be 'N' or 'iupac', got '{het_char}'")
    if lowcov not in ("n", "reference", "gnomad"):
        raise VflankError(f"--bam-lowcov must be n|reference|gnomad, got '{lowcov}'")
    return ConsensusPolicy(
        min_depth=min_depth, call_fract=call_fract, het_char=het_char,
        lowcov=lowcov, min_baseq=min_baseq, min_mapq=min_mapq,
    )


def load_bam_resolver(
    bam: Path | None, bam_map: Path | None
) -> tuple[Callable[[str], str | None] | None, int]:
    """Return ``(resolver, n_mapped)``.

    ``resolver(sample) -> bam_path | None``. With ``--bam`` the single BAM applies
    to every sample (single-sample runs); with ``--bam-map`` it is looked up by
    ``Tumor_Sample_Barcode``. Returns ``(None, 0)`` when no BAM is given.
    """
    if bam is not None and bam_map is not None:
        raise VflankError("Use either --bam or --bam-map, not both.")

    if bam is not None:
        if not bam.exists():
            raise VflankError(f"--bam not found: {bam}")
        path = str(bam)
        return (lambda _sample: path), -1  # -1 marks "single BAM, all samples"

    if bam_map is not None:
        if not bam_map.exists():
            raise VflankError(f"--bam-map not found: {bam_map}")
        mapping: dict[str, str] = {}
        for line in bam_map.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].strip() and parts[1].strip():
                mapping[parts[0].strip()] = parts[1].strip()
        if not mapping:
            raise VflankError(f"--bam-map produced no entries: {bam_map}")
        return (lambda sample: mapping.get(sample)), len(mapping)

    return None, 0
