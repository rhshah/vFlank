"""Shared CLI presentation helpers."""

from __future__ import annotations

from rich.panel import Panel
from rich.table import Table

from ..logging import console

# Help-panel group names for `--help`, shared across CLIs so the grouping is
# consistent (Typer's rich_help_panel). The most-used options (input, genome
# build, flank, output) stay in the default panel; these group the rest.
PANEL_REFERENCE = "Reference source"
PANEL_MASKING = "SNP masking (gnomAD)"
PANEL_BAM = "Patient consensus (BAM, modes C/D)"
PANEL_FILTER = "Sample selection"
PANEL_MAF_COLS = "MAF column names"
PANEL_SV_COLS = "Breakpoint column names"
PANEL_ADVANCED = "Advanced"


def echo_parameters(params: dict) -> None:
    """Print every effective run parameter as a panel before processing."""
    table = Table(show_header=False, box=None, padding=(0, 2, 0, 0))
    table.add_column(style="bold cyan", justify="right")
    table.add_column(overflow="fold")
    for key, value in params.items():
        table.add_row(key, str(value))
    console.print(Panel(table, title="[bold]Run parameters[/bold]", expand=False))
