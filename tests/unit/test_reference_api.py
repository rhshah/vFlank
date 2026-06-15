import os

import pytest

from vflank.core.reference_api import (
    ReferenceApiSource,
    _TransientApiError,
    build_url,
    genome_for_build,
    parse_sequence_response,
    ucsc_contig,
)
from vflank.errors import ReferenceError


def test_genome_for_build():
    assert genome_for_build("hg19") == "hg19"
    assert genome_for_build("hg38") == "hg38"
    with pytest.raises(ReferenceError):
        genome_for_build("hg99")


def test_ucsc_contig_prefix_and_mito():
    assert ucsc_contig("17") == "chr17"
    assert ucsc_contig("X") == "chrX"
    assert ucsc_contig("MT") == "chrM"  # not chrMT
    assert ucsc_contig("M") == "chrM"


def test_build_url_contains_genome_and_coords():
    url = build_url("hg19", "chr17", 100, 110)
    # UCSC coords are 0-based half-open — same as ours, no translation.
    assert "genome=hg19" in url and "chrom=chr17" in url
    assert "start=100" in url and "end=110" in url


# --- pure parser ---

def test_parse_returns_dna():
    assert parse_sequence_response({"dna": "ACGTacgt", "chrom": "chr17"}) == "ACGTacgt"


def test_parse_error_payload_raises():
    with pytest.raises(ReferenceError, match="boom"):
        parse_sequence_response({"error": "boom"})


def test_parse_missing_dna_raises():
    with pytest.raises(ReferenceError, match="missing a 'dna'"):
        parse_sequence_response({"chrom": "chr17", "start": 100})


def test_parse_short_sequence_is_not_an_error():
    # truncation at a contig end is a real condition the caller reports, not raised here.
    assert parse_sequence_response({"dna": "AC"}) == "AC"


# --- source with injected transport (no network) ---

class FakeTransport:
    def __init__(self, bodies):
        self.bodies = list(bodies)
        self.urls: list[str] = []

    def __call__(self, url, timeout):
        self.urls.append(url)
        b = self.bodies.pop(0)
        if isinstance(b, Exception):
            raise b
        return b


def _source(transport, **kw):
    return ReferenceApiSource(
        "hg19", transport=transport, sleep_fn=lambda _s: None, clock=lambda: 0.0, **kw,
    )


def test_fetch_coords_and_parse():
    t = FakeTransport([{"dna": "ACGTACGTAC"}])
    src = _source(t)
    assert src.fetch("17", 100, 110) == "ACGTACGTAC"
    assert "chrom=chr17" in t.urls[0] and "start=100" in t.urls[0] and "end=110" in t.urls[0]


def test_cache_avoids_duplicate_requests():
    t = FakeTransport([{"dna": "ACGT"}])
    src = _source(t)
    src.fetch("17", 100, 104)
    src.fetch("17", 100, 104)  # identical window -> cached
    assert src.request_count == 1
    assert len(t.urls) == 1


def test_empty_window_makes_no_request():
    t = FakeTransport([])
    src = _source(t)
    assert src.fetch("17", 50, 50) == ""
    assert src.request_count == 0


def test_throttle_sleeps_between_distinct_requests():
    waits = []
    src = ReferenceApiSource(
        "hg19",
        transport=FakeTransport([{"dna": "A"}, {"dna": "C"}]),
        sleep_fn=lambda s: waits.append(s), clock=lambda: 0.0, min_interval=1.0,
    )
    src.fetch("17", 100, 101)
    src.fetch("17", 200, 201)  # second distinct window
    assert waits and waits[0] == 1.0  # throttled to the rate limit


def test_retry_then_success():
    t = FakeTransport([_TransientApiError("timeout"), {"dna": "ACGT"}])
    src = _source(t)
    assert src.fetch("17", 100, 104) == "ACGT"
    assert len(t.urls) == 2


def test_transient_exhausted_raises_reference_error():
    t = FakeTransport([_TransientApiError("x")] * 10)
    src = _source(t, max_retries=2)
    with pytest.raises(ReferenceError, match="unreachable"):
        src.fetch("17", 100, 104)


def test_api_error_payload_raises():
    t = FakeTransport([{"error": "no such chromosome"}])
    src = _source(t)
    with pytest.raises(ReferenceError, match="no such chromosome"):
        src.fetch("17", 100, 104)


def test_check_build_trusts_declared():
    assert _source(FakeTransport([])).check_build("hg19") is None


@pytest.mark.skipif(
    not os.environ.get("VFLANK_LIVE_API"),
    reason="hits the live UCSC API; set VFLANK_LIVE_API=1 to run",
)
def test_live_api_tp53_window():
    src = ReferenceApiSource("hg19")
    seq = src.fetch("17", 7579400, 7579410)  # 10 bp window
    assert len(seq) == 10
    assert set(seq.upper()) <= set("ACGTN")
    src.close()
