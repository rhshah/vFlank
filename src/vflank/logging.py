"""Centralised Rich-backed logging and console setup.

A single stderr ``Console`` is shared across the package so that data written to
stdout (or a file) stays clean while human-facing status goes to stderr.
"""

from __future__ import annotations

import logging

from rich.console import Console
from rich.logging import RichHandler

# Status/console output always goes to stderr; stdout is reserved for data.
console = Console(stderr=True)

_LOGGER_NAME = "vflank"


def setup_logging(verbosity: int = 0, *, show_tracebacks: bool = False) -> logging.Logger:
    """Configure and return the ``vflank`` logger.

    Parameters
    ----------
    verbosity:
        ``-1`` = WARNING (quiet), ``0`` = INFO (default), ``>=1`` = DEBUG (verbose).
    show_tracebacks:
        When True, Rich renders rich tracebacks with locals on uncaught errors.
    """
    level = {
        -1: logging.WARNING,
        0: logging.INFO,
    }.get(verbosity, logging.DEBUG)

    handler = RichHandler(
        console=console,
        rich_tracebacks=show_tracebacks,
        tracebacks_show_locals=show_tracebacks,
        show_time=False,
        show_path=False,
        markup=True,
    )

    logger = logging.getLogger(_LOGGER_NAME)
    logger.handlers.clear()
    logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    return logger


def get_logger() -> logging.Logger:
    """Return the package logger (configured lazily at INFO if untouched)."""
    logger = logging.getLogger(_LOGGER_NAME)
    if not logger.handlers:
        setup_logging(0)
    return logger
