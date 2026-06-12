"""Population allele-frequency masking source via the gnomAD GraphQL API.

An alternative to :class:`~vflank.core.popfreq.GnomadStore` that needs no local
VCF download — it queries https://gnomad.broadinstitute.org/api per flank region.
Exposes the same duck-typed ``get_positions`` interface, so it drops in behind
``ReferenceFlankSource`` unchanged.

Trade-offs (see docs/research/gnomad-api.md): no download and both builds, but
rate-limited to ~10 requests/IP/60s and not reproducible — best for small
cohorts; prefer the VCF source for bulk/HPC/reproducible runs.

The parsing kernel (:func:`parse_api_variants`) is pure and unit-testable; HTTP
and timing are injected so tests run offline.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from collections.abc import Callable, Iterable
from time import monotonic, sleep

from ..errors import PopFreqError
from ..logging import get_logger
from .popfreq import kinds_for

log = get_logger()

API_URL = "https://gnomad.broadinstitute.org/api"

# vflank build -> (GraphQL reference_genome enum, dataset enum).
# gnomAD v4 is GRCh38-only; GRCh37 lives in v2.1.1.
DATASET_FOR_BUILD: dict[str, tuple[str, str]] = {
    "hg19": ("GRCh37", "gnomad_r2_1"),
    "hg38": ("GRCh38", "gnomad_r4"),
}


def dataset_for_build(genome_build: str) -> tuple[str, str]:
    try:
        return DATASET_FOR_BUILD[genome_build]
    except KeyError:
        raise PopFreqError(
            f"No gnomAD API dataset for build {genome_build!r} (expected hg19/hg38)."
        ) from None


def build_query(
    chrom: str, start: int, stop: int, reference_genome: str, dataset: str, kinds: Iterable[str]
) -> str:
    """GraphQL query for SNP AFs in a region. ``start``/``stop`` are 1-based inclusive."""
    seq = " ".join(f"{k} {{ af populations {{ id ac an }} }}" for k in kinds)
    return (
        f'{{ region(chrom: "{chrom}", start: {start}, stop: {stop}, '
        f"reference_genome: {reference_genome}) "
        f"{{ variants(dataset: {dataset}) {{ pos ref alt {seq} }} }} }}"
    )


def _variant_af(variant: dict, kinds: Iterable[str]) -> float:
    """Max AF for a variant across the chosen kinds: overall af and per-pop ac/an."""
    best = 0.0
    for kind in kinds:
        seq = variant.get(kind)
        if not seq:
            continue
        af = seq.get("af")
        if af is not None:
            best = max(best, af)
        for pop in seq.get("populations") or []:
            an = pop.get("an")
            if an:
                best = max(best, pop["ac"] / an)
    return best


def parse_api_variants(
    variants: Iterable[dict], kinds: Iterable[str], af_threshold: float
) -> list[int]:
    """1-based positions of SNPs whose max AF (over kinds) >= threshold.

    Pure: ``variants`` is the GraphQL ``region.variants`` list. SNPs only.
    """
    kinds = tuple(kinds)
    positions: list[int] = []
    for v in variants:
        if len(v.get("ref", "")) != 1 or len(v.get("alt", "")) != 1:
            continue
        if _variant_af(v, kinds) >= af_threshold:
            positions.append(int(v["pos"]))
    return positions


class _TransientApiError(Exception):
    """Retryable failure (network, timeout, 429, 5xx)."""


def _http_transport(url: str, query: str, timeout: float) -> dict:
    """POST a GraphQL query and return parsed JSON. Classifies failures."""
    data = json.dumps({"query": query}).encode()
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json"}
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (https only)
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        if exc.code == 429 or 500 <= exc.code < 600:
            raise _TransientApiError(f"HTTP {exc.code}") from exc
        raise PopFreqError(f"gnomAD API HTTP {exc.code}: {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise _TransientApiError(f"network error: {exc.reason}") from exc
    except json.JSONDecodeError as exc:
        raise PopFreqError(f"gnomAD API returned non-JSON: {exc}") from exc


class GnomadApiSource:
    """Masking source backed by the public gnomAD GraphQL API.

    Region responses are cached (so the two flank queries of identical variants
    reuse one request), requests are throttled to respect the rate limit, and
    transient failures are retried with backoff before raising ``PopFreqError``.
    """

    def __init__(
        self,
        genome_build: str,
        pop_data: str = "genome",
        *,
        url: str = API_URL,
        timeout: float = 30.0,
        min_interval: float = 6.0,   # ~10 requests / 60 s
        max_retries: int = 3,
        transport: Callable[[str, str, float], dict] | None = None,
        sleep_fn: Callable[[float], None] | None = None,
        clock: Callable[[], float] = monotonic,
    ) -> None:
        self.reference_genome, self.dataset = dataset_for_build(genome_build)
        self.kinds = kinds_for(pop_data)  # raises ValueError on bad input
        self.url = url
        self.timeout = timeout
        self.min_interval = min_interval
        self.max_retries = max_retries
        # Resolve defaults at init time (not as bound defaults) so tests can
        # monkeypatch the module-level transport/sleep.
        self._transport = transport if transport is not None else _http_transport
        self._sleep = sleep_fn if sleep_fn is not None else sleep
        self._clock = clock
        self._cache: dict[tuple[str, int, int], list[dict]] = {}
        self._last_call: float | None = None
        self.request_count = 0  # for monitoring / large-cohort warnings

    def _throttle(self) -> None:
        if self._last_call is not None:
            wait = self.min_interval - (self._clock() - self._last_call)
            if wait > 0:
                self._sleep(wait)
        self._last_call = self._clock()

    def _query_region(self, chrom: str, start: int, stop: int) -> list[dict]:
        query = build_query(chrom, start, stop, self.reference_genome, self.dataset, self.kinds)
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                self.request_count += 1
                body = self._transport(self.url, query, self.timeout)
            except _TransientApiError as exc:
                if attempt < self.max_retries:
                    self._sleep(2.0 * (attempt + 1))  # linear backoff
                    log.debug("gnomAD API transient error (%s), retrying", exc)
                    continue
                raise PopFreqError(
                    f"gnomAD API unreachable after {self.max_retries + 1} attempts: {exc}"
                ) from exc

            errors = body.get("errors")
            if errors:
                msg = errors[0].get("message", "unknown error")
                if "rate" in msg.lower() and attempt < self.max_retries:
                    self._sleep(2.0 * (attempt + 1))
                    log.debug("gnomAD API rate-limited, retrying")
                    continue
                raise PopFreqError(f"gnomAD API error: {msg}")

            variants = (((body.get("data") or {}).get("region") or {}).get("variants")) or []
            log.debug("gnomAD API %s:%d-%d -> %d variants", chrom, start, stop, len(variants))
            return variants
        return []  # unreachable; satisfies type checker

    def get_positions(
        self, bare: str, start_0based: int, end_0based: int, af_threshold: float
    ) -> list[int]:
        """1-based positions of common SNPs in ``[start, end)`` for this chrom."""
        if end_0based <= start_0based:
            return []
        start, stop = start_0based + 1, end_0based  # 0-based half-open -> 1-based inclusive
        key = (bare, start, stop)
        if key not in self._cache:
            self._cache[key] = self._query_region(bare, start, stop)
        return parse_api_variants(self._cache[key], self.kinds, af_threshold)

    def close(self) -> None:
        """No resources to release; present for interface symmetry with GnomadStore."""
