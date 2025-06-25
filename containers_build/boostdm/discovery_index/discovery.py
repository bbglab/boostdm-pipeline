import json
import tqdm
import glob
import os
import click
import random

import pandas as pd
import numpy as np

from boostdm.oncotree import Oncotree


def good_turing_discovery(freq_vector):
    """
    It computes the coverage, i.e. the total probability mass of the observed mutations,
    using the method by Good-Turing.
    
    Reference: Esty, W. W. A Normal Limit Law for a Nonparametric Estimator of the Coverage of a Random Sample. 
    The Annals of Statistics 11, 905-912 (1983).
    """

    f_one = np.sum(freq_vector == 1)
    n = np.sum(freq_vector)
    return 1 - (f_one / n)


@click.command()
@click.option('--evaluation-path', type=str)
@click.option('--mutations', type=click.Path(exists=True), required=True)
@click.option('--output', type=str)
def cli(evaluation_path, mutations, output):

    tree = Oncotree()

    df_discovery_index = {
        'gene': [], 
        'ttype': [], 
        'n_muts': [], 
        'n_unique_muts': [],
        'discovery_index': []
        }

    df_mutations = pd.read_csv(mutations, sep='\t')

    # filter SNVs
    df_mutations = df_mutations.loc[(df_mutations['REF'].isin(list('ACGT'))) & (df_mutations['ALT'].isin(list('ACGT')))]

    gene_ttype = {}
    file_iterable = list(glob.glob(os.path.join(evaluation_path, '*/*.eval.pickle.gz')))

    for fn in file_iterable:
        gene = os.path.basename(fn).split('.')[0]
        ttype = os.path.basename(os.path.dirname(fn))
        gene_ttype[ttype] = gene_ttype.get(ttype, []) + [gene]

    for ttype in tqdm.tqdm(gene_ttype):
        cohorts = tree.get_cohorts(ttype)
        df = df_mutations[df_mutations['COHORT'].isin(cohorts)]
        for gene in gene_ttype[ttype]:
            df_ungrouped = df[df['SYMBOL'] == gene]
            df_ungrouped['TUMOR_TYPE'] = ttype
            df_grouped = df_ungrouped.groupby(by=['MUTATION', 'SYMBOL', 'TUMOR_TYPE']).agg({'SAMPLES': 'sum'}).reset_index()
            freq_vector = df_grouped['SAMPLES'].values
            discovery_index = good_turing_discovery(freq_vector)
            # fill out the table
            df_discovery_index['gene'].append(gene)
            df_discovery_index['ttype'].append(ttype)
            df_discovery_index['n_muts'].append(np.sum(freq_vector))
            df_discovery_index['n_unique_muts'].append(len(freq_vector))
            df_discovery_index['discovery_index'].append(discovery_index)
        
    df_discovery_index = pd.DataFrame(df_discovery_index)
    df_discovery_index.to_csv(output, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    cli()
