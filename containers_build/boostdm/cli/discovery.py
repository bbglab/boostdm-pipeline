"""Cli commands for discovery."""
import click
from boostdm.discovery_index.muts import cli as prepare_dataset
from boostdm.discovery_index.preprocess_variants import cli as collect_metadata
from boostdm.discovery_index.samples import cli as discover_samples
from boostdm.discovery_index.discovery import cli as run

@click.group(name="discover")
def discovery_group():
    """Discover-related commands."""

discovery_group.add_command(prepare_dataset, name="prepare-dataset")
discovery_group.add_command(collect_metadata, name="collect-metadata")
discovery_group.add_command(discover_samples, name="collect-samples")
discovery_group.add_command(run, name="run")
