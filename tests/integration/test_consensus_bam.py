"""Synthetic-BAM tests for the samtools-consensus engine (Phase 1)."""

from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from vflank.core.consensus import (  # noqa: E402
    ConsensusPolicy,
    run_samtools_consensus,
    window_depth,
)

BASE = "ACGTG"
READ_LEN = 30
REF_START = 10  # 0-based; reads cover 0-based [10, 40)


def _template():
    return [BASE[i % len(BASE)] for i in range(READ_LEN)]


def _make_bam(tmp_path, *, n_reads=20, het_idx=None, het_alt="T"):
    """A BAM of identical reads, optionally heterozygous (50/50) at one index."""
    bam = tmp_path / "t.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 60}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(n_reads):
            seq = _template()
            if het_idx is not None and i >= n_reads // 2:
                seq[het_idx] = het_alt if seq[het_idx] != het_alt else "A"
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = REF_START
            a.mapping_quality = 60
            a.cigartuples = [(0, READ_LEN)]
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            out.write(a)
    pysam.index(str(bam))
    return str(bam)


def _genomic_to_read_idx(pos_1based):
    return pos_1based - REF_START - 1  # template index for a genomic position


def test_consensus_homozygous_matches_reads(tmp_path):
    bam = _make_bam(tmp_path)
    seq = run_samtools_consensus(bam, "chr1", 22, 32, ConsensusPolicy())
    assert len(seq) == 11
    # genomic pos 30 -> template base
    assert seq[30 - 22] == BASE[_genomic_to_read_idx(30) % len(BASE)]


def test_heterozygous_becomes_N(tmp_path):
    # het at genomic 25 (read idx 14: ref 'G' vs alt 'T') -> default consensus N
    bam = _make_bam(tmp_path, het_idx=14, het_alt="T")
    seq = run_samtools_consensus(bam, "chr1", 22, 32, ConsensusPolicy())
    assert seq[25 - 22] == "N"


def test_heterozygous_iupac(tmp_path):
    bam = _make_bam(tmp_path, het_idx=14, het_alt="T")
    seq = run_samtools_consensus(bam, "chr1", 22, 32, ConsensusPolicy(het_char="iupac"))
    assert seq[25 - 22] in "RYSWKM"  # an ambiguity code, not N


def test_window_depth(tmp_path):
    bam = _make_bam(tmp_path, n_reads=20)
    # 0-based [21, 32): all covered by 20 reads
    assert window_depth(bam, "chr1", 21, 32, ConsensusPolicy()) == [20] * 11
    # a region outside the read span -> zero depth
    assert window_depth(bam, "chr1", 45, 50, ConsensusPolicy()) == [0] * 5


def test_consensus_flank_source(tmp_path):
    from vflank.core.consensus import BamConsensusSource, ConsensusFlankSource
    from vflank.core.variant import Variant
    from vflank.io.reference import ReferenceFasta

    ref_seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{ref_seq}\n")
    pysam.faidx(str(fasta))

    rs, rl, nr = 25, 40, 24          # 24 reads cover 1-based 26..65
    het_g0 = 44                      # genomic 0-based 44 -> het (ref 'A' vs 'G')
    homalt_g0 = 54                   # genomic 0-based 54 -> hom-ALT (ref 'G' -> 'T')
    bam = tmp_path / "s.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(nr):
            seq = list(ref_seq[rs:rs + rl])
            seq[homalt_g0 - rs] = "T"             # all reads: hom-ALT
            if i >= nr // 2:
                seq[het_g0 - rs] = "G"            # half reads: het alt
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = rs
            a.mapping_quality = 60
            a.cigartuples = [(0, rl)]
            a.query_qualities = pysam.qualitystring_to_array("I" * rl)
            out.write(a)
    pysam.index(str(bam))

    ref = ReferenceFasta(str(fasta))
    src = ConsensusFlankSource(ref, BamConsensusSource(str(bam)), gnomad=None, flank=10)
    fr = src.fetch(Variant(chrom="1", start=50, end=50, ref="A", alt="T"))
    # left window 1-based 40..49; het at genomic 45 -> idx 5 -> N
    assert fr.masked_left[5] == "N"
    # right window 1-based 51..60; hom-ALT at genomic 55 -> idx 4 -> patient 'T' (ref 'G')
    assert fr.masked_right[4] == "T" and fr.right[4] == "G"
    # everything else equals the reference flank (raw)
    assert fr.masked_left[:5] + fr.masked_left[6:] == fr.left[:5] + fr.left[6:]
    ref.close()


def test_small_run_with_bam_consensus(tmp_path):
    from typer.testing import CliRunner

    from vflank.cli.app import app

    ref_seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{ref_seq}\n")
    pysam.faidx(str(fasta))

    rs, rl, nr = 25, 40, 24
    het_g0 = 44   # genomic 0-based -> het in the left flank
    bam = tmp_path / "S1.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(nr):
            seq = list(ref_seq[rs:rs + rl])
            if i >= nr // 2:
                seq[het_g0 - rs] = "G"
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = rs
            a.mapping_quality = 60
            a.cigartuples = [(0, rl)]
            a.query_qualities = pysam.qualitystring_to_array("I" * rl)
            out.write(a)
    pysam.index(str(bam))

    header_cols = "\t".join([
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ])
    maf = tmp_path / "v.maf"
    maf.write_text(header_cols + "\n" + "\t".join(["GENE", "1", "50", "50", "A", "T", "S1"]) + "\n")
    bammap = tmp_path / "map.tsv"
    bammap.write_text(f"S1\t{bam}\n")
    out_fa = tmp_path / "out.fasta"

    result = CliRunner().invoke(app, [
        "small", "run", str(maf), "-r", str(fasta), "-g", "hg38",
        "--bam-map", str(bammap), "-f", "10", "--output", str(out_fa),
    ])
    assert result.exit_code == 0, result.output
    lines = out_fa.read_text().splitlines()
    assert len(lines) == 4
    assert lines[0].startswith(">S1__GENE")            # sample back in the header
    assert lines[2].startswith(">Masked__S1__GENE")
    masked_left = lines[3].split("[", 1)[0]
    assert "N" in masked_left                          # het in left flank -> N in consensus
    assert "BAM consensus" in result.output


def test_consensus_handles_deletion_reads_reference_length(tmp_path):
    # Reads with a 3 bp deletion in the left flank must NOT change the consensus
    # length (regression: samtools indel output broke window alignment).
    from vflank.core.consensus import BamConsensusSource, ConsensusFlankSource
    from vflank.core.variant import Variant
    from vflank.io.reference import ReferenceFasta

    ref_seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{ref_seq}\n")
    pysam.faidx(str(fasta))

    rs, a, d, b, nr = 25, 15, 3, 22, 24   # M15 D3 M22 -> ref span 40, covers [25,65)
    query = ref_seq[rs:rs + a] + ref_seq[rs + a + d:rs + a + d + b]
    bam = tmp_path / "del.bam"
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        for i in range(nr):
            seg = pysam.AlignedSegment()
            seg.query_name = f"r{i}"
            seg.query_sequence = query
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = rs
            seg.mapping_quality = 60
            seg.cigartuples = [(0, a), (2, d), (0, b)]
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(query))
            out.write(seg)
    pysam.index(str(bam))

    ref = ReferenceFasta(str(fasta))
    src = ConsensusFlankSource(ref, BamConsensusSource(str(bam)), gnomad=None, flank=10)
    fr = src.fetch(Variant(chrom="1", start=50, end=50, ref="A", alt="T"))  # must not raise
    assert len(fr.masked_left) == 10 and len(fr.masked_right) == 10
    assert fr.total == 20  # reference-aligned throughout
    ref.close()


def test_consensus_flags_insertion_site(tmp_path):
    # Reads carrying a 4 bp insertion in the right flank: the consensus stays
    # reference-length (insertion dropped) but the anchor position is flagged as
    # N and reported via ``inserted`` — not lost silently.
    from vflank.core.consensus import (
        BamConsensusSource,
        ConsensusPolicy,
        insertion_sites,
    )
    from vflank.io.reference import ReferenceFasta

    ref_seq = "".join("ACGT"[i % 4] for i in range(100))
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">1\n{ref_seq}\n")
    pysam.faidx(str(fasta))

    rs, a, ins, b, nr = 49, 6, 4, 20, 24   # M6 I4 M20 -> ref span 26, covers [49,75)
    anchor0 = rs + a - 1                    # 0-based ref base before the insertion
    query = ref_seq[rs:rs + a] + "TTTT" + ref_seq[rs + a:rs + a + b]
    bam = tmp_path / "ins.bam"
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        for i in range(nr):
            seg = pysam.AlignedSegment()
            seg.query_name = f"r{i}"
            seg.query_sequence = query
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = rs
            seg.mapping_quality = 60
            seg.cigartuples = [(0, a), (1, ins), (0, b)]
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(query))
            out.write(seg)
    pysam.index(str(bam))

    policy = ConsensusPolicy(min_depth=10)
    assert insertion_sites(str(bam), "1", 49, 75, policy) == {anchor0}

    # Right flank of a variant at 50 is ref[50, 60) 0-based -> includes anchor0=54.
    ref = ReferenceFasta(str(fasta))
    bsrc = BamConsensusSource(str(bam), policy)
    seq, _covered, n_ins = bsrc.consensus("1", 50, 60, ref_seq[50:60], set())
    assert n_ins == 1
    assert seq[anchor0 - 50] == "N"
    assert len(seq) == 10  # reference-length preserved
    ref.close()
