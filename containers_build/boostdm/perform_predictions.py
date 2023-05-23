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

from boostdm.globals import COLUMNS_TRAINING, COLUMNS_OUTPUT, COLUMNS_SHAP, SYSTEMATIC_BIAS
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


def weighted_consensus(model_obj, eval_obj, use_weights=True):
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

    def func(x):

        prod = 1
        for i, model in enumerate(models):
            feature_names = model.model.get_booster().feature_names
            p = model.predict_proba(x[feature_names])[:, 1]
            if use_weights:
                prod *= (p / (1 - p)) ** weights[i]  # weights based on logloss
            else:
                prod *= (p / (1 - p)) ** (1 / len(weights))  # balanced weighting
        prod = prod ** SYSTEMATIC_BIAS  # systematic bias correction
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

        consensus_func = weighted_consensus(model, evaluation)
        mutations.loc[df.index, ['boostDM_score']] = consensus_func(df[COLUMNS_TRAINING])

        x_data = df[COLUMNS_TRAINING]
        shap_bootstrap = []
        for model in model['models']:
            explainer = shap.TreeExplainer(model.model)
            shap_bootstrap.append(explainer.shap_values(x_data, check_additivity=False))
        shap_values = np.mean(shap_bootstrap, axis=0)

        mutations.loc[df.index, COLUMNS_SHAP] = shap_values

    mutations.loc[:, 'boostDM_class'] = mutations['boostDM_score'].apply(lambda x: x >= 0.5)
    return mutations


def predict(mutations, gene, ttype, model_selection_dict, models_folder, evaluations_folder, high_quality_only=True):

    # model selection

    selected_model_ttype, _ = model_selection_dict.get((ttype, gene), (None, gene))

    # prepare mutations table

    for c in COLUMNS_SHAP:
        mutations[c] = np.nan
    mutations['boostDM_score'] = np.nan

    # select ttype model

    if high_quality_only and (selected_model_ttype is None):
        raise Exception("There is no gene-specific, high-quality model")
    elif (not high_quality_only) and (selected_model_ttype is None):
        selected_model_ttype = 'CANCER'

    # build consensus predictor

    path_model = os.path.join(models_folder, f'{selected_model_ttype}', f'{gene}.models.pickle.gz')
    with gzip.open(path_model, 'rb') as f:
        model = pickle.load(f)

    path_eval = os.path.join(evaluations_folder, f'{selected_model_ttype}', f'{gene}.eval.pickle.gz')
    with gzip.open(path_eval, 'rb') as g:
        evaluation = pickle.load(g)

    consensus_func = weighted_consensus(model, evaluation)

    # predict on mutations

    mutations.loc[mutations.index, ['boostDM_score']] = consensus_func(mutations[COLUMNS_TRAINING])

    # SHAP explanations

    x_data = mutations[COLUMNS_TRAINING]
    shap_bootstrap = []
    for model in model['models']:
        explainer = shap.TreeExplainer(model.model)
        shap_bootstrap.append(explainer.shap_values(x_data, check_additivity=False))
    shap_values = np.mean(shap_bootstrap, axis=0)

    mutations.loc[mutations.index, COLUMNS_SHAP] = shap_values
    mutations.loc[:, 'boostDM_class'] = mutations['boostDM_score'].apply(lambda x: x >= 0.5)

    mutations['selected_model_ttype'] = ttype

    return mutations


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
            model_selection_dict = pickle.load(g)
        df = predict(df_mutations,
                     gene=gene,
                     ttype=tumor_type,
                     model_selection_dict=model_selection_dict,
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
