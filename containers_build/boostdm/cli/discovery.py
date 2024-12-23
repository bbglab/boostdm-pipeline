"""Cli commands for discovery."""
import click
from boostdm.discovery_index.muts import cli as discover_muts
from boostdm.discovery_index.preprocess_variants import cli as collect
from boostdm.discovery_index.samples import cli as discover_samples
from boostdm.discovery_index.discovery import cli as discover

@click.group(name="discover")
def discovery_group():
    """Annotation-related commands."""
    pass

discovery_group.add_command(discover_muts, name="mutations")
discovery_group.add_command(collect, name="collect")
discovery_group.add_command(discover_samples, name="samples")
discovery_group.add_command(discover, name="run")
