"""Cli commands for annotations."""
import click
from boostdm.cvdata.cohort import cli as cohorts
from boostdm.cvdata.meta import cli as meta

@click.group(name="annotations")
def cvdata_group():
    """Perform cross-validation split operations."""

cvdata_group.add_command(cohorts, name="cohort")
cvdata_group.add_command(meta, name="meta")
