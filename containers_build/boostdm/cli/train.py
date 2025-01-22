"""Cli commands for training."""
import click
from boostdm.features.group import cli as group_features
from boostdm.annotations.cohort import cli as prepare_dataset
from boostdm.training import cli as model_train
from boostdm.evaluation.auto import cli as evaluate_model
from boostdm.cli._crossvalidation import cvdata_group

@click.group(name="train")
def cli():
    """BoostDM: Training-related commands."""

cli.add_command(group_features, name="group-features")
cli.add_command(prepare_dataset, name="prepare-dataset")
cli.add_command(cvdata_group, name="cv-splits")
cli.add_command(model_train, name="model-train")
cli.add_command(evaluate_model, name="evaluate-model")
