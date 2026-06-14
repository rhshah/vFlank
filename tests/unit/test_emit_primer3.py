from vflank.io.emit_primer3 import (
    junction_record,
    n_runs,
    resolve_template,
    small_variant_record,
    write_primer3,
)


def test_resolve_template_fills_N_from_reference():
    # masked/consensus call wins; reference base fills where the call is N.
    assert resolve_template("AAACC", "AGNCC") == "AGACC"


def test_n_runs_coalesces_and_offsets():
    assert n_runs("AANNCAN") == [(2, 2), (6, 1)]
    assert n_runs("ANA", offset=10) == [(11, 1)]
    assert n_runs("ACGT") == []


def test_small_variant_record_target_and_excluded():
    rec = small_variant_record(
        "TP53__x", left="AAACC", right="GGGTT",
        masked_left="AANCC", masked_right="GGGTN", ref="C",
    )
    assert rec.template == "AAACC" + "C" + "GGGTT"   # N filled from reference
    assert rec.target == (5, 1)                       # variant 'C' after the left flank
    assert rec.excluded == [(2, 1), (10, 1)]          # masked sites, right shifted by 6


def test_junction_record_target_straddles_junction():
    rec = junction_record("FUS__x", "AAAATTTT", "AAANTTTT", junction_index=4)
    assert rec.template == "AAAATTTT"
    assert rec.target == (3, 2)         # last base of partner1 + first of partner2
    assert rec.excluded == [(3, 1)]


def test_write_primer3_boulder_io_structure(tmp_path):
    recs = [
        small_variant_record("v1", "AAA", "TTT", "ANA", "TTT", "C"),
        small_variant_record("v2", "GGG", "CCC", "GGG", "CCC", "A"),
    ]
    out = tmp_path / "p3.txt"
    write_primer3(out, recs)
    text = out.read_text()

    # Global settings appear once, in the first record only.
    assert text.count("PRIMER_PICK_INTERNAL_OLIGO=1") == 1
    assert "PRIMER_PRODUCT_SIZE_RANGE=60-200" in text
    # Two records, each terminated by a lone '='.
    assert text.count("\n=\n") == 2
    assert "SEQUENCE_ID=v1" in text and "SEQUENCE_ID=v2" in text
    assert "SEQUENCE_TARGET=3,1" in text                 # variant after the 3 bp left flank
    assert "SEQUENCE_EXCLUDED_REGION=1,1" in text        # the masked N in v1's left flank
    # v2 has no masked positions -> no excluded-region tag in its block.
    v2_block = text.split("SEQUENCE_ID=v2")[1]
    assert "SEQUENCE_EXCLUDED_REGION" not in v2_block
