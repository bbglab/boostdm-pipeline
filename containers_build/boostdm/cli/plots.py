"""Cli commands for plots."""
import click
from boostdm.output_plots.blueprint import cli as blueprint
from boostdm.output_plots.clustered_blueprint import cli as clustered_blueprint
from boostdm.output_plots.discovery_plot import cli as discovery_plot

@click.group(name="plot")
def cli():
    """BoostDM: Plotting-related commands."""
    pass

cli.add_command(blueprint, name="blueprint")
cli.add_command(clustered_blueprint, name="clustered_blueprint")
cli.add_command(discovery_plot, name="discovery")
