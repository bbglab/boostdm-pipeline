import click
from boostdm.features.group import cli as group_features
from boostdm.annotations.cohort import cli as prepare_dataset
from boostdm.training import cli as model_train
from boostdm.evaluation.auto import cli as evaluate_model
from boostdm.cli._crossvalidation import cvdata_group

@click.group(name="train")
def train_group():
    """Training-related commands."""

train_group.add_command(group_features, name="group-features")
train_group.add_command(prepare_dataset, name="prepare-dataset")
train_group.add_command(cvdata_group, name="cv-splits")
train_group.add_command(model_train, name="model-train")
train_group.add_command(evaluate_model, name="evaluate-model")
