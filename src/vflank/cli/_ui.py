"""Shared CLI presentation helpers."""

from __future__ import annotations

from rich.panel import Panel
from rich.table import Table

from ..logging import console


def echo_parameters(params: dict) -> None:
    """Print every effective run parameter as a panel before processing."""
    table = Table(show_header=False, box=None, padding=(0, 2, 0, 0))
    table.add_column(style="bold cyan", justify="right")
    table.add_column(overflow="fold")
    for key, value in params.items():
        table.add_row(key, str(value))
    console.print(Panel(table, title="[bold]Run parameters[/bold]", expand=False))
