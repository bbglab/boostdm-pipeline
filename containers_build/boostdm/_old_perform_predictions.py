"""Script to conduct prediction with gradient boosting models"""

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

from boostdm.globals import COLUMNS_TRAINING, COLUMNS_OUTPUT, COLUMNS_SHAP
from boostdm.evaluation.data import Hierarchy

warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
warnings.filterwarnings(module='sklearn*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='pandas*', action='ignore', category=RuntimeWarning)


class ExtendedHierarchy(Hierarchy):

    def __init__(self, models):
        super().__init__()
        self._models = models

    def get_model(self, gene, ttype):
        return self._models.get((ttype, gene), (None, gene))

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

    # weights for consensus combination
    weights = 1 / np.array(logloss)
    weights = weights / np.sum(weights)

    # systematic bias
    bias = 2.3

    def func(x):

        prod = 1
        for i, model in enumerate(models):
            feature_names = model.model.get_booster().feature_names
            p = model.predict_proba(x[feature_names])[:, 1]
            if use_weights:
                prod *= (p / (1 - p)) ** weights[i]  # weights based on logloss
            else:
                prod *= (p / (1 - p)) ** (1 / len(weights))  # balanced weighting
        prod = prod ** bias  # systematic bias correction
        s = (prod / (1 + prod))

        return s

    return func


def _predict(mutations, models_folder, evaluations_folder, high_quality_only=True):

    for c in COLUMNS_SHAP:
        mutations[c] = np.nan
    mutations['boostDM_score'] = np.nan

    for name, df in mutations.groupby(['selected_model_gene', 'selected_model_ttype']):
        gene, ttype = name

        if high_quality_only:
            if ttype is None:
                raise Exception("There is no gene-specific, high-quality model")
        else:
            ttype = 'CANCER'

        path_model = os.path.join(models_folder, f'{ttype}', f'{gene}.models.pickle.gz')
        with gzip.open(path_model, 'rb') as f:
            model = pickle.load(f)

        path_eval = os.path.join(evaluations_folder, f'{ttype}', f'{gene}.eval.pickle.gz')
        with gzip.open(path_eval, 'rb') as g:
            evaluation = pickle.load(g)

        consensus_func = bootstrap_voting(model, evaluation)
        mutations.loc[df.index, ['boostDM_score']] = consensus_func(df[COLUMNS_TRAINING])

        x_data = df[COLUMNS_TRAINING]
        shap_bootstrap = []
        for model in model['models']:
            explainer = shap.TreeExplainer(model.model)
            shap_bootstrap.append(explainer.shap_values(x_data))
        shap_values = np.mean(shap_bootstrap, axis=0)

        mutations.loc[df.index, COLUMNS_SHAP] = shap_values

    mutations.loc[:, 'boostDM_class'] = mutations['boostDM_score'].apply(lambda x: x >= 0.5)
    return mutations


def predict(mutations, ttype, models, models_folder, evaluations_folder, high_quality_only=True):

    hierarchy = ExtendedHierarchy(models)
    _get_model = functools.partial(hierarchy.get_model_for_gene, ttype=ttype)
    d = dict(map(_get_model, mutations['gene'].unique()))
    mutations['selected_model_ttype'], mutations['selected_model_gene'] = zip(*mutations['gene'].map(d))

    return _predict(mutations, models_folder, evaluations_folder, high_quality_only=high_quality_only)


def predict_with_model(mutations, model, models_folder, evaluations_folder):

    ttype, gene = model
    mutations['selected_model_ttype'] = ttype
    mutations['selected_model_gene'] = gene

    return _predict(mutations, models_folder, evaluations_folder)


@click.command()
@click.option('--muts', type=click.Path(exists=True), help="File with the annotated mutations", required=True)
@click.option('--gene', type=str, help="gene symbol", required=True)
@click.option('--tumor-type', type=str, help="Tumor type of the input mutations", default=None)
@click.option('--models-folder', type=click.Path(exists=True), help="Path to the folder where the models are", required=True)
@click.option('--evaluations-folder', type=click.Path(exists=True), help="Path to the folder where the evaluations are", required=True)
@click.option('--model-selection', type=str, help="Either the path to the model selection dict, or the gene to be used", required=True)
@click.option('--output-file', type=click.Path(), help="File with the annotated mutations", required=True)
@click.option('--high-quality-only', is_flag=True, show_default=True, default=False, help="Prediction only if matching models are high-quality")
def cli(muts, gene, tumor_type, models_folder, evaluations_folder, model_selection, output_file, high_quality_only):
    
    df = pd.read_csv(muts, sep='\t')
    df_mutations = df[df['gene'] == gene]


    if os.path.exists(model_selection):
        with gzip.open(model_selection, 'rb') as g:
            models = pickle.load(g)
        df = predict(df_mutations,
                     ttype=tumor_type,
                     models=models,
                     models_folder=models_folder,
                     evaluations_folder=evaluations_folder,
                     high_quality_only=high_quality_only)
    else:
        df = predict_with_model(df_mutations, (tumor_type, model_selection),
                                models_folder=models_folder,
                                evaluations_folder=evaluations_folder)

    output_cols = COLUMNS_OUTPUT + COLUMNS_SHAP

    df[output_cols].to_csv(output_file,
                           index=False,
                           compression="gzip",
                           sep="\t")


if __name__ == '__main__':
    cli()
