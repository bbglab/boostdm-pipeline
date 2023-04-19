"""Script to conduct training of gradient boosting models"""

# Imports
import functools
import gzip
import os
import warnings

import click
import dill as pickle
import numpy as np
import pandas as pd
import shap

warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
warnings.filterwarnings(module='sklearn*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='pandas*', action='ignore', category=RuntimeWarning)

from boostdm.globals import COLUMNS_TRAINING
from boostdm.evaluation.data import Hierarchy as ModelHierarchy

COLUMNS_OUTPUT = ['gene', 'ENSEMBL_TRANSCRIPT', 'ENSEMBL_GENE', 'chr', 'pos', 'alt', 'aachange'] + \
                 COLUMNS_TRAINING + \
                 ['selected_model_ttype', 'selected_model_gene', 'boostDM_score', 'boostDM_class']
COLUMNS_SHAP = [f'shap_{x}' for x in COLUMNS_TRAINING]


class Hierarchy(ModelHierarchy):

    def __init__(self, models):
        super().__init__()
        self._models = models

    # TODO: define a final rule for model selection

    """
    def get_model(self, gene, ttype):
        if (ttype, gene) in self._models:
            return self._models[(ttype, gene)]
        else:
            for t in self.climb(ttype, gene):
                if t in self._models:
                    self._models[(ttype, gene)] = t
                    return self._models[t]
    """

    # TODO: replace get_model by its its final form

    @staticmethod
    def get_model(gene, ttype):
        return ttype, gene

    def get_model_for_gene(self, gene, ttype):
        return gene, self.get_model(gene, ttype)


def bootstrap_voting(model_obj, eval_obj, use_weights=True):
    """
    Args:
        model_obj: instance of model for a specific ttype-gene
        eval_obj: evaluation instance for the same specific ttype-gene
        use_weights
    Return:
        function: that given a vector of features, it returns a prediction label and score
                  based on individual predictions and voting among all models and thresholds
        warning:  at this point probabilities may need a final recalibration step
    Satopaa, V. A. et al. Combining multiple probability predictions using a simple logit model.
    International Journal of Forecasting 30, 344–356 (2014).
    """

    models = model_obj['models']
    logloss = eval_obj['logloss']

    # estimate overall systematic bias
    weights = 1 / np.array(logloss)
    weights = weights / np.sum(weights)
    bias = 2.3  # estimate of systematic bias inferred from BRCA TCGA data

    def func(x):

        prod = 1
        for i, model in enumerate(models):
            feature_names = model.model.get_booster().feature_names
            p = model.predict_proba(x[feature_names])[:, 1]
            if use_weights:
                prod *= (p / (1 - p)) ** weights[i]  # weighting based on logloss
            else:
                prod *= (p / (1 - p)) ** (1 / len(weights))  # otherwise, balanced weighting
        prod = prod ** bias  # correction for systematic bias
        s = (prod / (1 + prod))

        return s

    return func


def _predict_group(df, models_folder, evaluations_folder):
    gene, ttype = df.name

    path_model = os.path.join(models_folder, f'{ttype}', f'{gene}.models.pickle.gz')
    with gzip.open(path_model, 'rb') as f:
        model = pickle.load(f)

    path_eval = os.path.join(evaluations_folder, f'{ttype}', f'{gene}.eval.pickle.gz')
    with gzip.open(path_eval, 'rb') as g:
        evaluation = pickle.load(g)

    return bootstrap_voting(model, evaluation)(df[COLUMNS_TRAINING])


def _predict(mutations, models_folder, evaluations_folder):

    # Prepare space for the columns
    for c in COLUMNS_SHAP:
        mutations[c] = np.nan
    mutations['boostDM_score'] = np.nan

    for name, df in mutations.groupby(['selected_model_gene', 'selected_model_ttype']):
        gene, ttype = name

        path_model = os.path.join(models_folder, f'{ttype}', f'{gene}.models.pickle.gz')
        with gzip.open(path_model, 'rb') as f:
            model = pickle.load(f)

        path_eval = os.path.join(evaluations_folder, f'{ttype}', f'{gene}.eval.pickle.gz')
        with gzip.open(path_eval, 'rb') as g:
            evaluation = pickle.load(g)

        mutations.loc[df.index, ['boostDM_score']] = bootstrap_voting(model, evaluation)(df[COLUMNS_TRAINING])

        # Attribution
        x_data = df[COLUMNS_TRAINING]

        shap_bootstrap = []
        for model in model['models']:
            explainer = shap.TreeExplainer(model.model)
            shap_bootstrap.append(explainer.shap_values(x_data))
        shap_values = np.mean(shap_bootstrap, axis=0)

        mutations.loc[df.index, COLUMNS_SHAP] = shap_values

    mutations.loc[:, 'boostDM_class'] = mutations['boostDM_score'].apply(lambda x: x >= 0.5)
    return mutations


def predict(mutations, ttype, models, models_folder, evaluations_folder):

    hierarchy = Hierarchy(models)

    get_model = functools.partial(hierarchy.get_model_for_gene, ttype=ttype)

    d = dict(map(get_model, mutations['gene'].unique()))
    mutations['selected_model_ttype'], mutations['selected_model_gene'] = zip(*mutations['gene'].map(d))

    return _predict(mutations, models_folder, evaluations_folder)


def predict_with_model(mutations, model, models_folder, evaluations_folder):
    ttype, gene = model
    mutations['selected_model_ttype'] = ttype
    mutations['selected_model_gene'] = gene

    return _predict(mutations, models_folder, evaluations_folder)


def _attribute_group(df, models_folder):
    """
    d: dict of the form {ttype: [genes]}
    df: dataframe with features
    """
    gene, ttype = df.name

    path_model = os.path.join(models_folder, f'{ttype}', f'{gene}.models.pickle.gz')
    with gzip.open(path_model, 'rb') as f:
        model = pickle.load(f)

    x_data = df[COLUMNS_TRAINING]

    shap_bootstrap = []
    for model in model['models']:
        explainer = shap.TreeExplainer(model.model)
        shap_bootstrap.append(explainer.shap_values(x_data))
    shap_values = np.mean(shap_bootstrap, axis=0)
    return pd.DataFrame(shap_values, columns=COLUMNS_SHAP)


@click.command()
@click.option('--muts', type=click.Path(exists=True), help="File with the annotated mutations", required=True)
@click.option('--output-file', type=click.Path(), help="File with the annotated mutations", required=True)
@click.option('--tumor-type', type=str, help="Tumor type of the input mutations", default=None)
@click.option('--models-folder', type=click.Path(exists=True), help="Path to the folder where the models are", required=True)
@click.option('--evaluations-folder', type=click.Path(exists=True), help="Path to the folder where the evaluations are", required=True)
@click.option('--model-selection', type=str, help="Either the path to the model selection dict, or the gene to be used", required=True)
def cli(muts, output_file, tumor_type, model_selection, models_folder, evaluations_folder):
    df_mutations = pd.read_csv(muts, sep='\t')

    if os.path.exists(model_selection):
        # if there is a model selection file
        with gzip.open(model_selection, 'rb') as g:
            models = pickle.load(g)
        df = predict(df_mutations, ttype=tumor_type, models=models,
                     models_folder=models_folder, evaluations_folder=evaluations_folder)
    else:
        df = predict_with_model(df_mutations, (tumor_type, model_selection),
                                models_folder=models_folder, evaluations_folder=evaluations_folder)

    output_cols = COLUMNS_OUTPUT + COLUMNS_SHAP

    df[output_cols].to_csv(
        output_file,
        index=False,
        compression="gzip",
        sep="\t"
    )


if __name__ == '__main__':
    cli()
