import os

import pytest

from vflank.core.popfreq_api import (
    GnomadApiSource,
    _TransientApiError,
    build_query,
    dataset_for_build,
    parse_api_variants,
)
from vflank.errors import PopFreqError


def test_dataset_for_build():
    assert dataset_for_build("hg19") == ("GRCh37", "gnomad_r2_1")
    assert dataset_for_build("hg38") == ("GRCh38", "gnomad_r4")
    with pytest.raises(PopFreqError):
        dataset_for_build("hg99")


def test_build_query_contains_fields_and_coords():
    q = build_query("17", 101, 110, "GRCh37", "gnomad_r2_1", ("genome", "exome"))
    assert 'chrom: "17"' in q and "start: 101" in q and "stop: 110" in q
    assert "reference_genome: GRCh37" in q and "dataset: gnomad_r2_1" in q
    assert "genome {" in q and "exome {" in q


# --- pure parser ---

def _variants():
    return [
        {"pos": 105, "ref": "A", "alt": "G",
         "genome": {"af": 0.2, "populations": [{"id": "nfe", "ac": 10, "an": 50}]},
         "exome": {"af": 0.0, "populations": []}},
        {"pos": 107, "ref": "AT", "alt": "A",                 # indel -> excluded
         "genome": {"af": 0.9, "populations": []}, "exome": None},
        {"pos": 109, "ref": "C", "alt": "T",                  # rare -> excluded
         "genome": {"af": 0.0001, "populations": [{"id": "afr", "ac": 1, "an": 2000}]},
         "exome": None},
        {"pos": 111, "ref": "G", "alt": "A",                  # common only via exome
         "genome": {"af": 0.0, "populations": []},
         "exome": {"af": 0.0, "populations": [{"id": "nfe", "ac": 600, "an": 1000}]}},
    ]


def test_parse_genome_only():
    assert parse_api_variants(_variants(), ("genome",), 0.001) == [105]


def test_parse_both_picks_up_exome_pop():
    # pos 111 is common only in the exome population (ac/an=0.6) -> union includes it.
    assert parse_api_variants(_variants(), ("genome", "exome"), 0.001) == [105, 111]


# --- source with injected transport (no network) ---

class FakeTransport:
    def __init__(self, bodies):
        self.bodies = list(bodies)
        self.queries: list[str] = []

    def __call__(self, url, query, timeout):
        self.queries.append(query)
        b = self.bodies.pop(0)
        if isinstance(b, Exception):
            raise b
        return b


def _ok_body(variants):
    return {"data": {"region": {"variants": variants}}}


def _source(transport, **kw):
    return GnomadApiSource(
        "hg19", kw.pop("pop_data", "genome"),
        transport=transport, sleep_fn=lambda _s: None, clock=lambda: 0.0, **kw,
    )


def test_get_positions_coords_and_parse():
    t = FakeTransport([_ok_body(_variants())])
    src = _source(t)
    # 0-based half-open [100, 110) -> 1-based inclusive start=101 stop=110
    assert src.get_positions("17", 100, 110, 0.001) == [105]
    assert "start: 101" in t.queries[0] and "stop: 110" in t.queries[0]


def test_region_cache_avoids_duplicate_requests():
    t = FakeTransport([_ok_body(_variants())])
    src = _source(t)
    src.get_positions("17", 100, 110, 0.001)
    src.get_positions("17", 100, 110, 0.001)  # identical region -> cached
    assert src.request_count == 1
    assert len(t.queries) == 1


def test_empty_window_makes_no_request():
    t = FakeTransport([])
    src = _source(t)
    assert src.get_positions("17", 50, 50, 0.001) == []
    assert src.request_count == 0


def test_throttle_sleeps_between_distinct_requests():
    waits = []
    src = GnomadApiSource(
        "hg19", "genome",
        transport=FakeTransport([_ok_body([]), _ok_body([])]),
        sleep_fn=lambda s: waits.append(s), clock=lambda: 0.0, min_interval=6.0,
    )
    src.get_positions("17", 100, 110, 0.001)
    src.get_positions("17", 200, 210, 0.001)  # second distinct region
    assert waits and waits[0] == 6.0  # throttled to the rate limit


def test_retry_then_success():
    t = FakeTransport([_TransientApiError("timeout"), _ok_body(_variants())])
    src = _source(t)
    assert src.get_positions("17", 100, 110, 0.001) == [105]
    assert len(t.queries) == 2


def test_graphql_error_raises():
    t = FakeTransport([{"errors": [{"message": "boom"}]}])
    src = _source(t)
    with pytest.raises(PopFreqError, match="boom"):
        src.get_positions("17", 100, 110, 0.001)


def test_transient_exhausted_raises_popfreq():
    t = FakeTransport([_TransientApiError("x")] * 10)
    src = _source(t, max_retries=2)
    with pytest.raises(PopFreqError, match="unreachable"):
        src.get_positions("17", 100, 110, 0.001)


@pytest.mark.skipif(
    not os.environ.get("VFLANK_LIVE_API"),
    reason="hits the live gnomAD API; set VFLANK_LIVE_API=1 to run",
)
def test_live_api_tp53_p72r():
    src = GnomadApiSource("hg19", "genome")
    pos = src.get_positions("17", 7579440, 7579540, 0.01)
    assert 7579472 in pos  # TP53 P72R common SNP
    src.close()
