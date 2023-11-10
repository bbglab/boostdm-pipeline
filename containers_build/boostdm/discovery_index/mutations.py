import os
import json
import random
import tqdm
import glob

import pandas as pd
import numpy as np

from boostdm.oncotree import Oncotree
from boostdm.config import PathConfig

# annotmuts data from which we collect
base_dirs = [PathConfig.intogen_data]
cohorts_path = PathConfig.cohorts_path


def oncotree_children(ttype):
    """Generator of cohorts belonging to the same ttype type"""

    tree = Oncotree()
    cohorts_list = tree.get_cohorts(ttype)
    cohorts, projects = list(zip(*cohorts_list))
    return cohorts


def get_mutations(ttype):
    cohorts = oncotree_children(ttype)
    mutations_table = []
    for cohort in cohorts:
        for dir_ in base_dirs:
            path = os.path.join(dir_, f'dndscv/{cohort}.annotmuts.gz')
            if not os.path.exists(path):
                continue
            df = pd.read_csv(path, sep='\t')
            mutations_table.append(df)
    df = pd.concat(mutations_table, axis=0)
    df['chr'] = df['chr'].astype(str)
    df['pos'] = df['pos'].astype(int)
    return df


# variant sample information paths
samples_intogen = os.path.join(intogen, 'vep', 'variants.json')


def get_sample_info():
    df_stats = pd.read_csv(cohorts_path, sep="\t")
    df_stats.rename(columns={"CANCER_TYPE": "cancer_type"}, inplace=True)

    # samples info
    intogen_samples = json.load(open(samples_intogen, 'r'))

    samples_info = {}
    for d in [intogen_samples]:
        for cohort in d:
            try:
                ttype = df_stats[df_stats['COHORT'] == cohort]['cancer_type'].values[0]
            except:
                continue

            if not ttype in samples_info:
                samples_info[ttype] = []

            list_hypermutators = []
            if 'hypermutators' in d[cohort]:
                list_hypermutators += list(d[cohort]['hypermutators']['hypermutators'])
            samples_info[ttype] += list(d[cohort]['donors'].keys()) + list_hypermutators
    return samples_info


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
        l = []
        for _ in range(iterations):
            selected_samples = random.sample(samples_info, int(value))
            muts = df_observed[(df_observed['sampleID'].isin(selected_samples)) & (
                        df_observed['impact'] != 'Synonymous')][['chr', 'pos', 'mut']].drop_duplicates().shape[0]
            l.append(muts)
        unique_counts.append(l)
    return grid, unique_counts


def get_number_mutations(evaluation_folder):
    n_mutations = {}
    for fn in tqdm.tqdm(glob.glob(os.path.join(evaluation_folder, '*'))):
        ttype = os.path.basename(fn)
        df_observed = get_mutations(ttype)
        group = df_observed.groupby(by='gene')['sampleID'].count()
        for gn in list(glob.glob(os.path.join(evaluation_folder, f'{ttype}/*.eval.pickle.gz'))):
            gene = os.path.basename(gn).split('.')[0]
            if gene in ['LoF', 'Act', 'pan']:
                continue
            try:
                n_mutations[(gene, ttype)] = group[gene]
            except:
                print(gene, ttype)
    return n_mutations
