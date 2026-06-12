"""vflank root CLI.

    vflank small run ...        # extract + mask flanks from a MAF
    vflank small inspect ...    # preview MAF columns
    vflank small list-vcf ...   # verify gnomAD directory coverage
    vflank version
"""

from __future__ import annotations

import typer

from .. import __version__
from ..logging import setup_logging
from . import small

app = typer.Typer(
    name="vflank",
    help="Variant-aware flanking-sequence extraction and masking for ddPCR assay design.",
    add_completion=False,
    no_args_is_help=True,
)
app.add_typer(small.app, name="small", help="Small-variant (SNP/indel) flank extraction.")


@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable DEBUG logging."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Only show warnings and errors."),
    debug: bool = typer.Option(False, "--debug", help="DEBUG logging + rich tracebacks."),
):
    """Global options applied before any subcommand."""
    verbosity = 1 if (verbose or debug) else (-1 if quiet else 0)
    setup_logging(verbosity, show_tracebacks=debug)


@app.command()
def version():
    """Print the vflank version."""
    typer.echo(__version__)


if __name__ == "__main__":
    app()
