"""Cli commands for training."""

import click
from src.features.group import cli as group_features
from src.annotations.cohort import cli as prepare_dataset
from src.training import cli as model_train
from src.evaluation.auto import cli as evaluate_model
from src.cli._crossvalidation import cvdata_group


@click.group(name="train")
def cli():
    """BoostDM: Training-related commands."""


cli.add_command(group_features, name="group-features")
cli.add_command(prepare_dataset, name="prepare-dataset")
cli.add_command(cvdata_group, name="cv-splits")
cli.add_command(model_train, name="model-train")
cli.add_command(evaluate_model, name="evaluate-model")
