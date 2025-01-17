"""Cli commands for annotations."""
import click
from boostdm.annotations.gene import cli as prepare_dataset
from boostdm.perform_predictions import cli as saturation

@click.group(name="predict")
def predict_group():
    """Prediction-related commands."""

predict_group.add_command(prepare_dataset, name="prepare-dataset")
predict_group.add_command(saturation, name="saturation")
