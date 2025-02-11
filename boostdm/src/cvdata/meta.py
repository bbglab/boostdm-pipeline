"""Prepare CV-data for subsequent training"""

import gzip
import os
import pickle
import warnings
from collections import defaultdict
from multiprocessing import Pool

import click
import pandas as pd

from src.cvdata.utils import sort_filter, vertical_join
from src.oncotree import Oncotree

warnings.filterwarnings(module="sklearn*", action="ignore", category=DeprecationWarning)
warnings.filterwarnings(module="sklearn*", action="ignore", category=RuntimeWarning)
warnings.filterwarnings(module="pandas*", action="ignore", category=RuntimeWarning)


def generate(arg):
    """
    Aggregates cross-validation (CV) results for gene-related data from multiple input files.

    This function reads multiple gzip-compressed pickle files, each containing a dictionary
    where keys are gene names, and values are lists of tuples of pandas DataFrames. Each tuple
    represents the results of a cross-validation split. The function merges corresponding CV
    splits across files, preserving the structure while avoiding empty DataFrames.
    """
    _ttype, input_files, output_folder = arg

    dict_split_output = {}

    for file in input_files:
        with gzip.open(file, "rb") as f:
            d_cvobj = pd.read_pickle(f)
        for gene in d_cvobj:
            if d_cvobj[gene] is None:
                continue
            if gene in dict_split_output:
                dict_split_output[gene] = vertical_join(
                    dict_split_output[gene], d_cvobj[gene]
                )
            else:
                dict_split_output[gene] = d_cvobj[gene]

        for gene, data in d_cvobj.items():
            if not data or all(
                df.empty for df_tuple in data for df in df_tuple
            ):  # Skip if all are empty
                continue

            if gene in dict_split_output:
                # Merge corresponding CV split DataFrames across input files
                dict_split_output[gene] = [
                    tuple(
                        pd.concat([l, r], sort=False, axis=0).reset_index(drop=True)
                        for l, r in zip(existing, new)
                    )
                    for existing, new in zip(dict_split_output[gene], data)
                ]
            else:
                dict_split_output[gene] = data

    for gene, cvdata in dict_split_output.items():
        output_file = os.path.join(output_folder, f"{gene}.cvdata.pickle.gz")
        # sort filter the aggregated cross-validation splits
        cvdata = [sort_filter(*arg) for arg in cvdata]

        with gzip.open(output_file, "wb") as g:
            pickle.dump(cvdata, g)

    return


@click.command()
@click.option("--input_path", type=str)
@click.option("--output_path", type=str)
@click.option("--cores", type=int, default=None)
def cli(input_path, output_path, cores):
    """
    Args:
        input_path: path of the tables
        output_path: path where to save the outputs
        cores: number of cores to use
    Returns:
        pickle dump for cvdata.{ttype}.{gene}.pickle.gz
    """

    # TODO oncotree - it should be used the bbglab/oncotree package
    oncotree = Oncotree()
    ttypes = oncotree.ttypes
    cohorts = defaultdict(list)
    for ttype in ttypes:
        for cohort in oncotree.get_cohorts(ttype):
            cohort_file = os.path.join(input_path, f"{cohort}.cvdata.pickle.gz")
            if os.path.exists(cohort_file):
                cohorts[ttype].append(cohort_file)

    execution_tuples = []

    for ttype, input_files in cohorts.items():
        output_folder = os.path.join(output_path, f"{ttype}")
        os.makedirs(output_folder, exist_ok=True)
        # include ttype to save the time that it will take to send the data back to the main process
        execution_tuples.append((ttype, input_files, output_folder))

    with Pool(cores) as p:
        for _ in p.imap(generate, execution_tuples):
            pass


if __name__ == "__main__":
    cli()
