"""import modules"""
import click
from boostdm.cli.train import train_group
from boostdm.cli.discovery import discovery_group
from boostdm.evaluation.data import cli as evaluate_model
from boostdm.cli.predict import predict_group
from boostdm.cli.plots import plots_group
from boostdm.cli.benchmarks import benchmarks_group

@click.group()
def cli():
    """BoostDM: A tool for identifying driver mutations in cancer."""

# Add all groups to the main CLI
cli.add_command(train_group)
cli.add_command(discovery_group)
cli.add_command(evaluate_model, name='evaluate')
cli.add_command(predict_group)
cli.add_command(plots_group)
cli.add_command(benchmarks_group)
