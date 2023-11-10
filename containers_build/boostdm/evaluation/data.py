import gzip
import os
import pickle
import glob
import json

import click
import numpy as np
import pandas as pd

from boostdm.oncotree import Oncotree
from boostdm.globals import DISCOVERY_TIERS, MUTATION_TIERS, FSCORE_THRESHOLD


def get_fscore(model_evaluation):
    return np.nanmedian(model_evaluation.get('fscore50', 0.5))


def get_precision(model_evaluation):
    return np.nanmedian(model_evaluation.get('precision', 0.5))


def get_recall(model_evaluation):
    return np.nanmedian(model_evaluation.get('recall', 0.5))


def meet_condition(fscore, discovery, n_muts):
    if fscore >= FSCORE_THRESHOLD:
        for discovery_thresh, n_muts_thresh in zip(DISCOVERY_TIERS, MUTATION_TIERS):
            if (discovery >= discovery_thresh) and (n_muts >= n_muts_thresh):
                return True
    return False


def evaluate(model_evaluations):

    res = {}

    for ttype, gene in model_evaluations.keys():
        
        fscore = model_evaluations[(ttype, gene)]['fscore50']
        discovery = model_evaluations[(ttype, gene)]['discovery']
        n_muts = model_evaluations[(ttype, gene)]['n_muts']

        if meet_condition(fscore, discovery, n_muts):
            res[str((ttype, gene))] = str((ttype, gene))
        else:
            res[str((ttype, gene))] = str((None, None))
    return res

   
@click.command()
@click.option('--eval_folder', 'eval_folder', type=click.Path(), help='input folder containing autoevaluation results')
@click.option('--discovery_path', 'discovery_path', type=click.Path(), help='file path to discovery output table')
@click.option('--output', 'output_file', type=click.Path(), help='output file')
def cli(eval_folder, discovery_path, output_file):
    models = {}
    df_discovery = pd.read_csv(discovery_path, sep='\t')
    discovery_dict = df_discovery.set_index(['gene', 'ttype']).to_dict()['discovery_index']
    n_muts_dict = df_discovery.set_index(['gene', 'ttype']).to_dict()['n_muts']

    for fn in glob.glob(os.path.join(eval_folder, '*/*.eval.pickle.gz')):
        gene = os.path.basename(fn).split('.')[0]
        ttype = os.path.basename(os.path.dirname(fn))
        with gzip.open(fn, 'rb') as fd:
            d = pickle.load(fd)
        fscore50 = get_fscore(d)
        discovery = discovery_dict.get((gene, ttype), 0)
        n_muts = n_muts_dict.get((gene, ttype), 0)

        models[(ttype, gene)] = {'fscore50': fscore50, 'discovery': discovery, 'n_muts': n_muts}
    res = evaluate(models)
    with open(output_file, 'wt') as f:
        json.dump(res, f)


if __name__ == '__main__':
    cli()
