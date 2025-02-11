"""Cli commands for benchmarking."""
import click
from src.benchmarks import prepare_vep_input 
from src.benchmarks.create_cv_tables import cli as create_cv_tables
from src.benchmarks.saturation_dbnsfp import cli as saturation_dbnsfp
from src.benchmarks import annotate_cv_tables
from src.benchmarks.precision_recall import cli as precision_recall

@click.group(name="benchmark")
def cli():
    """BoostDM: Benchmark-related commands."""

cli.add_command(prepare_vep_input, name="prepare_vep")
cli.add_command(create_cv_tables, name="create")
cli.add_command(annotate_cv_tables, name="annotate")
cli.add_command(saturation_dbnsfp, name="saturation")
cli.add_command(precision_recall, name="run")
