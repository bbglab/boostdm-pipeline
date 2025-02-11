"""Cli commands for prediction."""
import click
from src.annotations.gene import cli as prepare_dataset
from src.perform_predictions import cli as saturation

@click.group(name="predict")
def cli():
    """BoostDM: Prediction-related commands."""

cli.add_command(prepare_dataset, name="prepare-dataset")
cli.add_command(saturation, name="saturation")
