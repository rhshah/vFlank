"""Typed exception hierarchy for vflank.

Keeping a small hierarchy lets the CLI map domain errors to clean, user-facing
messages while library callers can catch ``VflankError`` to handle anything.
"""

from __future__ import annotations


class VflankError(Exception):
    """Base class for all vflank errors."""


class ReferenceError(VflankError):
    """Problem with the reference FASTA or its index."""


class MafError(VflankError):
    """Problem reading or validating the input MAF."""


class PopFreqError(VflankError):
    """Problem locating or reading population-frequency (gnomAD) data."""
