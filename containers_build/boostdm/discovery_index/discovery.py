import json
import tqdm
import glob
import os
import click
import random

import pandas as pd
import numpy as np
from scipy.optimize import minimize

from boostdm.oncotree import Oncotree


def get_downsampling_counts(samples_info, df_observed, iterations=10, n_grid=10):
    """
    returns:
    - grid: list of values, representing number of samples from which unique counts are drawn
    - unique_counts: list of unique count set sizes for various replicates; each set has n=iterations elements
      the list has same length as n_grid
    """

    n = len(samples_info)
    grid = np.linspace(0, n, n_grid)
    unique_counts = []
    for value in grid:
        muts_list = []
        for _ in range(iterations):
            selected_samples = random.sample(samples_info, int(value))
            muts = df_observed[
                (df_observed['sampleID'].isin(selected_samples)) & (df_observed['impact'] != 'Synonymous')
            ][['chr', 'pos', 'mut']].drop_duplicates().shape[0]
            muts_list.append(muts)
        unique_counts.append(muts_list)
    return grid, unique_counts


def curve_fit(func, x_data, y_data, seed_params=None, weights=None, bounds=None):
    """custom curve fit function with a vector of weights"""

    if weights is None:
        weights = np.ones_like(x_data)
    assert(len(weights) == len(x_data))

    def cost(params):
        predictions = np.array([func(x, *params) for x in x_data])
        return np.dot(weights, (predictions - y_data) ** 2)

    res = minimize(cost, seed_params, bounds=bounds, options={'disp': False})
    return cost(res.x), res.x


def master_func(x, m, p):
    """
    master function to be fit
    refer to the Supplementary Note for a justification
    """
    return m * (1 - ((m - p) / m) ** x)


def discovery_index(total_samples, *params):
    """how far best fitting curve plateaus from current level unique mutations discovered"""
    m = params[0]
    a = master_func(total_samples, *params) / m
    return min(a, 1)


def bootstrap_data(unique_counts, iterations=10, ngrid=7):
    y_data_bootstrap = []
    for _ in range(iterations):
        y_sample = []
        for i in range(ngrid):
            r = np.random.choice(unique_counts[i], size=1)
            y_sample.append(r[0])
        y_data_bootstrap += [y_sample]
    return y_data_bootstrap


def fitting_with_bootstrap(grid, unique_counts, iterations=10, ngrid=7):
    params_pool = []
    y_data_bootstrap = bootstrap_data(unique_counts, iterations=iterations, ngrid=ngrid)

    # setting weights
    n = 2 * np.log10(ngrid)
    weights = 1 / (1 + np.std(y_data_bootstrap, axis=0) ** n)

    # setting optimization constraints
    bounds = [(unique_counts[-1][-1] / 2, None), (0, None)]

    for i, y in enumerate(y_data_bootstrap):
        error, params = curve_fit(master_func, grid, y,
                                  seed_params=[y[-1], 0.05],
                                  weights=weights,
                                  bounds=bounds)
        params_pool += [list(params)]

    total_samples = grid[-1]
    sat_list = [discovery_index(total_samples, *params) for params in params_pool]

    return params_pool, sat_list


def discovery_index_with_bootstrap(samples, mutations, iterations=10, ngrid=10):
    grid, unique_counts = get_downsampling_counts(samples, mutations,
                                                  iterations=iterations, n_grid=ngrid)
    params, disc_ind = fitting_with_bootstrap(grid, unique_counts, iterations=iterations, ngrid=ngrid)
    return params, disc_ind, grid, unique_counts


def discovery_run(samples, mutations, iterations=100, ngrid=20):

    np.random.seed(123)
    random.seed(123)

    params_list, disc_ind, grid, unique_counts = discovery_index_with_bootstrap(samples, mutations, iterations, ngrid)
    median = np.nanmedian(disc_ind)
    interquartile_range = (np.nanquantile(disc_ind, 0.25), np.nanquantile(disc_ind, 0.75))
    return grid[-1], unique_counts[-1][-1], median, interquartile_range


@click.command()
@click.option('--evaluation-path', type=str)
@click.option('--mutations', type=click.Path(exists=True), required=True)
@click.option('--samples', type=click.Path(exists=True), required=True)
@click.option('--output', type=str)
def cli(evaluation_path, mutations, samples, output):

    tree = Oncotree()

    df_discovery_index = {'gene': [], 'ttype': [], 'n_muts': [], 'n_unique_muts': [],
                          'n_samples': [], 'discovery_index': [], 'discovery_high': [], 'discovery_low': []}

    df_mutations = pd.read_csv(mutations, sep='\t')
    with open(samples, 'r') as fd:
        samples_info = json.load(fd)

    gene_ttype_iterable = {}
    file_iterable = list(glob.glob(os.path.join(evaluation_path, '*/*.eval.pickle.gz')))

    for fn in file_iterable:
        gene = os.path.basename(fn).split('.')[0]
        ttype = os.path.basename(os.path.dirname(fn))
        gene_ttype_iterable[ttype] = gene_ttype_iterable.get(ttype, []) + [gene]

    for ttype in tqdm.tqdm(gene_ttype_iterable):
        if ttype not in samples_info:
            continue
        cohorts = tree.get_cohorts(ttype)
        df_observed_ttype = df_mutations[df_mutations['COHORT'].isin(cohorts)]
        for gene in gene_ttype_iterable[ttype]:
            df_observed = df_observed_ttype[df_observed_ttype['gene'] == gene]
            n_muts = df_observed['sampleID'].count()
            try:
                n_samples, n_unique, discovery, interquartile_range = discovery_run(samples_info[ttype], df_observed)
            except Exception:
                #TODO: implement a more consistent logging for failed cases?
                print(f'Discovery index calculation failed for gene {gene}, tumor-type {ttype}')
                continue
            df_discovery_index['gene'].append(gene)
            df_discovery_index['ttype'].append(ttype)
            df_discovery_index['n_muts'].append(n_muts)
            df_discovery_index['n_unique_muts'].append(n_unique)
            df_discovery_index['n_samples'].append(int(n_samples))
            df_discovery_index['discovery_index'].append(discovery)
            df_discovery_index['discovery_high'].append(interquartile_range[1])
            df_discovery_index['discovery_low'].append(interquartile_range[0])
        
    df_discovery_index = pd.DataFrame(df_discovery_index)
    # df_discovery_index = df_discovery_index[~df_discovery_index.isnull()]
    df_discovery_index.to_csv(output, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    cli()
