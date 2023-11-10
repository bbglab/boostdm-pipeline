"""Prepare CV-data for subsequent training"""

# TODO py file: rename as meta

import gzip
import pickle
import os
from collections import defaultdict
from multiprocessing import Pool

import click


import warnings
warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
warnings.filterwarnings(module='sklearn*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='pandas*', action='ignore', category=RuntimeWarning)

from boostdm.cvdata.utils import vertical_join, sort_filter
from boostdm.oncotree import Oncotree


def generate(arg):
    ttype, input_files, output_folder = arg

    dict_split_output = {}

    for file in input_files:
        with gzip.open(file, 'rb') as f:
            d_cvobj = pickle.load(f)
        for gene in d_cvobj:
            if d_cvobj[gene] is None:
                continue
            if gene in dict_split_output:
                dict_split_output[gene] = vertical_join(dict_split_output[gene], d_cvobj[gene])
            else:
                dict_split_output[gene] = d_cvobj[gene]

    for gene, cvdata in dict_split_output.items():
        output_file = os.path.join(output_folder, f'{gene}.cvdata.pickle.gz')
        # sort filter the aggregated cross-validation splits
        cvdata = [sort_filter(*arg) for arg in cvdata]

        with gzip.open(output_file, 'wb') as g:
            pickle.dump(cvdata, g)

    return


@click.command()
@click.option('--input_path', type=str)
@click.option('--output_path', type=str)
@click.option('--cores', type=int, default=None)
def cli(input_path, output_path, cores):
    """
    Args:
        input_path: path of the tables
        output_path: path where to save the outputs
        cores: number of cores to use
    Returns:
        pickle dump for cvdata.{ttype}.{gene}.pickle.gz
    """

    oncotree = Oncotree()
    ttypes = oncotree.ttypes
    cohorts = defaultdict(list)
    for ttype in ttypes:
        for cohort, _ in oncotree.get_cohorts(ttype):
            cohort_file = os.path.join(input_path, f'{cohort}.cvdata.pickle.gz')
            if os.path.exists(cohort_file):
                cohorts[ttype].append(cohort_file)

    if not os.path.exists(output_path):  #TODO: warning, temporary hack!
        os.makedirs(output_path, exist_ok=True)

        execution_tuples = []
        for ttype, input_files in cohorts.items():
            output_folder = os.path.join(output_path, f'{ttype}')
            os.makedirs(output_folder, exist_ok=True)
            # include ttype to save the time that it will take to send the data back to the main process
            execution_tuples.append((ttype, input_files, output_folder))

        with Pool(cores) as p:
            for _ in p.imap(generate, execution_tuples):
                pass


if __name__ == '__main__':
    cli()
