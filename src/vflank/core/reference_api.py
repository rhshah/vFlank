"""Reference-sequence source backed by the UCSC REST API.

An alternative to :class:`~vflank.io.reference.ReferenceFasta` that needs no
local FASTA download — it fetches each flank window from
``https://api.genome.ucsc.edu/getData/sequence``. Exposes the surface the CLIs
use on a reference (``fetch`` / ``check_build`` / ``has_chr`` / ``close``), so it
drops in behind :class:`~vflank.core.flanks.ReferenceFlankSource` unchanged.

Why UCSC (see docs/research/genome-api.md): its genome ids are literally
``hg19`` / ``hg38`` (our ``--genome-build`` values) and its coordinates are
0-based half-open — identical to pysam — so the flank math in ``flanks.py`` needs
no translation. Trade-offs: no download and both builds, but a ~1 req/s courtesy
rate limit and network required — best for the hosted single-variant / small
use; prefer a local FASTA for bulk/HPC/offline runs.

The parsing kernel (:func:`parse_sequence_response`) and URL builder
(:func:`build_url`) are pure and unit-testable; HTTP and timing are injected so
tests run offline.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from time import monotonic, sleep

from ..errors import ReferenceError
from ..logging import get_logger

log = get_logger()

API_URL = "https://api.genome.ucsc.edu/getData/sequence"

# vflank build -> UCSC `genome` id. UCSC names match our build strings exactly.
GENOME_FOR_BUILD: dict[str, str] = {"hg19": "hg19", "hg38": "hg38"}


def genome_for_build(genome_build: str) -> str:
    try:
        return GENOME_FOR_BUILD[genome_build]
    except KeyError:
        raise ReferenceError(
            f"No UCSC genome for build {genome_build!r} (expected hg19/hg38)."
        ) from None


def ucsc_contig(bare: str) -> str:
    """Bare chromosome -> UCSC contig name.

    UCSC hg19/hg38 use ``chr``-prefixed names, and the mitochondrion is ``chrM``
    (not ``chrMT``). We cannot cheaply probe the contig list over the API, so we
    apply UCSC's known convention rather than guessing.
    """
    if bare in ("MT", "M"):
        return "chrM"
    return f"chr{bare}"


def build_url(genome: str, chrom: str, start_0based: int, end_0based: int) -> str:
    """UCSC getData/sequence URL. UCSC coords are 0-based half-open — same as ours."""
    query = urllib.parse.urlencode(
        {"genome": genome, "chrom": chrom, "start": start_0based, "end": end_0based}
    )
    return f"{API_URL}?{query}"


def parse_sequence_response(body: dict) -> str:
    """Extract the sequence from a UCSC getData/sequence JSON body.

    Pure. Raises :class:`ReferenceError` on an API error payload or a missing/
    non-string ``dna`` field. A *shorter-than-requested* sequence is **not** an
    error here — like a pysam fetch off a contig end, truncation is a real
    condition the caller reports (see the flank-truncation check in the CLIs).
    """
    if isinstance(body.get("error"), str):
        raise ReferenceError(f"UCSC API error: {body['error']}")
    dna = body.get("dna")
    if not isinstance(dna, str):
        raise ReferenceError(
            f"UCSC API response missing a 'dna' string (got keys: {sorted(body)})"
        )
    return dna


class _TransientApiError(Exception):
    """Retryable failure (network, timeout, 429, 5xx)."""


def _http_transport(url: str, timeout: float) -> dict:
    """GET a UCSC API URL and return parsed JSON. Classifies failures."""
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (https only)
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        if exc.code == 429 or 500 <= exc.code < 600:
            raise _TransientApiError(f"HTTP {exc.code}") from exc
        raise ReferenceError(f"UCSC API HTTP {exc.code}: {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise _TransientApiError(f"network error: {exc.reason}") from exc
    except json.JSONDecodeError as exc:
        raise ReferenceError(f"UCSC API returned non-JSON: {exc}") from exc


class ReferenceApiSource:
    """Reference-sequence source backed by the UCSC getData/sequence API.

    Window responses are cached, requests are throttled to respect UCSC's
    courtesy limit (~1 req/s), and transient failures are retried with backoff
    before raising ``ReferenceError``. Drop-in for :class:`ReferenceFasta` across
    the surface the CLIs use.
    """

    # UCSC serves chr-prefixed contigs; surfaced for the CLI's status line only.
    has_chr = True

    def __init__(
        self,
        genome_build: str,
        *,
        url: str = API_URL,
        timeout: float = 30.0,
        min_interval: float = 1.0,   # ~1 request / s (UCSC courtesy limit)
        max_retries: int = 3,
        transport: Callable[[str, float], dict] | None = None,
        sleep_fn: Callable[[float], None] | None = None,
        clock: Callable[[], float] = monotonic,
    ) -> None:
        self.genome = genome_for_build(genome_build)  # raises on bad build
        self.genome_build = genome_build
        self.url = url
        self.timeout = timeout
        self.min_interval = min_interval
        self.max_retries = max_retries
        # Resolve defaults at init time (not as bound defaults) so tests can
        # monkeypatch the module-level transport/sleep.
        self._transport = transport if transport is not None else _http_transport
        self._sleep = sleep_fn if sleep_fn is not None else sleep
        self._clock = clock
        self._cache: dict[tuple[str, int, int], str] = {}
        self._last_call: float | None = None
        self.request_count = 0  # for monitoring

    def _throttle(self) -> None:
        if self._last_call is not None:
            wait = self.min_interval - (self._clock() - self._last_call)
            if wait > 0:
                self._sleep(wait)
        self._last_call = self._clock()

    def _request(self, contig: str, start_0based: int, end_0based: int) -> str:
        url = build_url(self.genome, contig, start_0based, end_0based)
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                self.request_count += 1
                body = self._transport(url, self.timeout)
            except _TransientApiError as exc:
                if attempt < self.max_retries:
                    self._sleep(2.0 * (attempt + 1))  # linear backoff
                    log.debug("UCSC API transient error (%s), retrying", exc)
                    continue
                raise ReferenceError(
                    f"UCSC API unreachable after {self.max_retries + 1} attempts: {exc}"
                ) from exc
            dna = parse_sequence_response(body)
            log.debug("UCSC API %s:%d-%d -> %d bp", contig, start_0based, end_0based, len(dna))
            return dna
        return ""  # unreachable; satisfies the type checker

    def fetch(self, bare: str, start_0based: int, end_0based: int) -> str:
        """Reference bases for ``[start, end)`` (0-based half-open), bare chromosome."""
        if end_0based <= start_0based:
            return ""
        contig = ucsc_contig(bare)
        key = (contig, start_0based, end_0based)
        if key not in self._cache:
            self._cache[key] = self._request(contig, start_0based, end_0based)
        return self._cache[key]

    def check_build(self, declared: str) -> str | None:
        """No local sequence to fingerprint; the API serves the requested build.

        Returns ``None`` (no mismatch warning is possible). We trust ``declared``
        and pass it through as the UCSC ``genome``; a wrong build still surfaces
        downstream as a UCSC error rather than silent wrong sequence.
        """
        log.debug("Reference build %r trusted (UCSC API serves the requested genome)", declared)
        return None

    def close(self) -> None:
        """No resources to release; present for interface symmetry with ReferenceFasta."""
