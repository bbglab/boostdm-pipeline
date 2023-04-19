"""Prepare CV-data for subsequent training"""
import gzip
import pickle

import click
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit

import warnings
warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
warnings.filterwarnings(module='sklearn*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='pandas*', action='ignore', category=RuntimeWarning)

from boostdm.cvdata.utils import sort_filter, vertical_join


def split_balanced(x_data, y_data, test_size=0.3):
    """Generate balanced train-test split"""

    one_index  = list(y_data[y_data == 1].index)
    zero_index = list(y_data[y_data == 0].index)

    # randomly select n_ones indices from zero_index
    zero_index = list(np.random.choice(zero_index, size=len(one_index), replace=False))

    x_data_sub = x_data.loc[one_index + zero_index, :]
    y_data_sub = y_data.loc[one_index + zero_index]

    # the random state should be fixed prior to this call
    x_train, x_test, y_train, y_test = train_test_split(x_data_sub, y_data_sub,
                                                        test_size=test_size)
    return x_train, x_test, y_train, y_test


def get_cv_sets_balanced(x_data, y_data, n, size):
    """Generate several balanced train-test sets"""

    for _ in range(n):
        x_train, x_test, y_train, y_test = split_balanced(x_data, y_data, test_size=size)
        yield x_train, x_test, y_train, y_test


def prepare(data, nsplits=10, test_size=0.2):

    data.rename(columns={'response': 'label'}, inplace=True)

    # keep 'pos' and 'chr' to run position-based filtering
    # avoid = ['cohort', 'gene', 'ref', 'alt', 'aachange', 'label', 'motif']
    avoid = ['cohort', 'gene', 'aachange', 'label', 'motif']  # include 'ref' 'alt'
    features = list(filter(lambda x: x not in avoid, data.columns))

    # regression datasets
    x_data = data[features]
    y_data = data['label']

    # cv_list output has tuples (x_train, x_test, y_train, y_test) as elements
    cv_list = get_cv_sets_balanced(x_data, y_data, nsplits, test_size)

    # filter test data to prevent site repetitions
    cv_list = [sort_filter(*arg) for arg in cv_list]

    return cv_list


def combine_genes(name, genes, d_output):
    combination_result = None

    for gene in genes:
        cvlist = d_output[gene]
        if combination_result is None:
            combination_result = cvlist
        else:
            combination_result = vertical_join(combination_result, cvlist)

    d_output[name] = combination_result


def load_mutations(input_path):
    data = pd.read_csv(input_path, sep='\t', low_memory=False)
    return data


def generate(mutations, random_state=None, bootstrap_splits=50, cv_fraction=0.3):
    # fix random seed (if provided)
    np.random.seed(random_state)

    d_output = {}

    for gene, data in mutations.groupby('gene'):
        cv_list = prepare(data.copy(), nsplits=bootstrap_splits, test_size=cv_fraction)
        d_output[gene] = cv_list

    # add meta-cohorts
    # genes_act = mutations[mutations['role_Act'] == 1]['gene'].unique()
    # combine_genes('Act', genes_act, d_output)
    # genes_lof = mutations[mutations['role_LoF'] == 1]['gene'].unique()
    # combine_genes('LoF', genes_lof, d_output)
    # genes_pan = mutations['gene'].unique()
    # combine_genes('pan', genes_pan, d_output)

    return d_output


@click.command()
@click.option('--input_path', type=str)
@click.option('--output_path', type=str)
@click.option('--splits', 'bootstrap_splits', default=50, type=int)
@click.option('--seed', 'random_state', default=None, type=int)
@click.option('--cv', 'cv_fraction', default=0.3, type=float)
def cli(input_path, output_path, random_state, bootstrap_splits, cv_fraction):

    mutations = load_mutations(input_path)

    d_output = generate(mutations, random_state=random_state,
                                  bootstrap_splits=bootstrap_splits, cv_fraction=cv_fraction)

    with gzip.open(output_path, 'wb') as f:
        pickle.dump(d_output, f)


if __name__ == '__main__':
    cli()
