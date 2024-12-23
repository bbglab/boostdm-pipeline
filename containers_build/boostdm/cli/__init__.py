"""import modules"""
import click
from boostdm.cli.annotations import annotations_group
from boostdm.cli.discovery import discovery_group
from boostdm.cli.plots import plots_group
from boostdm.cli.benchmarks import benchmarks_group

@click.group()
def cli():
    """BoostDM: A tool for identifying driver mutations in cancer."""
    pass

# Add all groups to the main CLI
cli.add_command(annotations_group)
cli.add_command(discovery_group)
cli.add_command(plots_group)
cli.add_command(benchmarks_group)
