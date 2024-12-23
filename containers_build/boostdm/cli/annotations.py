"""Cli commands for annotations."""
import click
from boostdm.annotations.cohort import cli as annotate_cohort
from boostdm.annotations.gene import cli as annotate_gene

@click.group(name="annotations")
def annotations_group():
    """Annotation-related commands."""
    pass

annotations_group.add_command(annotate_cohort, name="cohort")
annotations_group.add_command(annotate_gene, name="gene")


