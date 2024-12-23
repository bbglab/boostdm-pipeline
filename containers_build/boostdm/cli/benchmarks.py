"""Cli commands for benchmarking."""
import click
from boostdm.benchmarks import prepare_vep_input 
from boostdm.benchmarks.create_cv_tables import cli as create_cv_tables
from boostdm.benchmarks.saturation_dbnsfp import cli as saturation_dbnsfp
from boostdm.benchmarks import annotate_cv_tables
from boostdm.benchmarks.precision_recall import cli as precision_recall

@click.group(name="benchmark")
def benchmarks_group():
    """Annotation-related commands."""
    pass

benchmarks_group.add_command(prepare_vep_input, name="prepare_vep")
benchmarks_group.add_command(create_cv_tables, name="create")
benchmarks_group.add_command(annotate_cv_tables, name="annotate")
benchmarks_group.add_command(saturation_dbnsfp, name="saturation")
benchmarks_group.add_command(precision_recall, name="run")
