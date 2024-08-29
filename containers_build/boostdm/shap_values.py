"""
For each mutation produce its shapley-additive values (shap-values) arising from its selected model.

Do the shap-value representation by doing aggregation of shap-values by feature class:
- Cluster 1D
- Cluster 3D
- Domain enrichment
- Conservation
- Role
- PTMs
- Consequence type

For each model, take the pool of all the mutations for which the model is selected.
From this shap-value matrix we can infer an overall feature attribution for the model.

"""


import sys
#sys.path.append('./..')

import gzip
import pickle
import os
import warnings

import click
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
import tqdm


warnings.filterwarnings('ignore')

from boostdm.oncotree import Oncotree


"""
Constants
"""

oncotree = Oncotree()
all_cohorts = oncotree.get_cohorts('CANCER')


def config_params(font_size=7):
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


@click.group()
def cli():
    pass


def bootstrap_shap(model_object, x_data):
    """
    Given a model instance -- a suite of models -- and a collection
    of feature values, it generates the shapley value predictions.
    """

    shap_bootstrap = []
    for model in model_object['models']:
        explainer = shap.TreeExplainer(model.model)
        shap_bootstrap.append(explainer.shap_values(x_data, check_additivity=False))
    return np.mean(shap_bootstrap, axis=0)


def attribution_func(d, df, models_path):

    """
    d: dict of the form {ttype: [genes]}
    df: dataframe with features
    """

    ttype = next(iter(d.keys()))

    with gzip.open(os.path.join(models_path, f'{ttype}.models.pickle.gz'), 'rb') as f:
        m = pickle.load(f)

    df_pool = []
    for gene in d[ttype]:

        model_object = m[gene]
        features = list(model_object['models'][0].model.get_booster().feature_names)
        df_slice = df[(df['selected_model_ttype'] == ttype) & (df['selected_model_gene'] == gene)].copy()
        x_data = df_slice[features].copy()

        new_columns = list(map(lambda x: 'shap_' + x, features))
        for col in new_columns:
            df_slice[col] = None
        shapley_values = bootstrap_shap(model_object, x_data)
        df_slice.loc[:, new_columns] = shapley_values

        df_pool.append(df_slice)

    dg = pd.concat(df_pool, axis=0)

    return dg


def get_attribution(df, models_path):

    models_iterable = list(set(df.apply(lambda row: (row['selected_model_ttype'],
                                                     row['selected_model_gene']),
                                        axis=1)))
    d = {}
    for ttype, gene in models_iterable:
        d[ttype] = d.get(ttype, []) + [gene]
    iterable = [{k: v} for k, v in d.items()]
    dataframe_pool = []
    for d in tqdm.tqdm(iterable):
        dg = attribution_func(d, df, models_path)
        dataframe_pool.append(dg)
    res = pd.concat(dataframe_pool, axis=0)

    return res


@cli.command()
@click.option('--input', type=click.Path(), help='table of results after attribution run')
@click.option('--output', type=click.Path(), help='output table with new shap_columns')
def shapley_heatmap(input, output):

    df = pd.read_csv(input, sep='\t')
    df = df[[c for c in df.columns if c.startswith('shap_group_')]]
    params = {'vmin': 0, 'vmax': 2, 'col_cluster': False, 'row_cluster': False,
              'cmap': 'YlOrRd', 'xticklabels': True, 'yticklabels': False, 'figsize': (10, 10)}
    g = sns.clustermap(df, **params)
    g.ax_heatmap.set_xlabel('feature class')
    g.ax_heatmap.set_ylabel('sample')

    g.savefig(output, dpi=300, bbox_inches='tight', transparent=False)
    plt.show()


if __name__ == '__main__':
    cli()
