"""Cli commands for discovery."""
import click
from src.discovery_index.muts import cli as prepare_dataset
from src.discovery_index.preprocess_variants import cli as collect_metadata
from src.discovery_index.samples import cli as discover_samples
from src.discovery_index.discovery import cli as run

@click.group(name="discover")
def cli():
    """boostdm: Discover-related commands."""

cli.add_command(prepare_dataset, name="prepare-dataset")
cli.add_command(collect_metadata, name="collect-metadata")
cli.add_command(discover_samples, name="collect-samples")
cli.add_command(run, name="run")
