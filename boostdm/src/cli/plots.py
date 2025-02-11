"""Cli commands for plots."""
import click
from src.output_plots.blueprint import cli as blueprint
from src.output_plots.clustered_blueprint import cli as clustered_blueprint
from src.output_plots.discovery_plot import cli as discovery_plot

@click.group(name="plot")
def cli():
    """boostdm: Plotting-related commands."""
    pass

cli.add_command(blueprint, name="blueprint")
cli.add_command(clustered_blueprint, name="clustered_blueprint")
cli.add_command(discovery_plot, name="discovery")
