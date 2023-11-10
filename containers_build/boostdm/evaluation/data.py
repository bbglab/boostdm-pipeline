import glob
import gzip
import os
import pickle

import click
import numpy as np
import pandas as pd

from boostdm.oncotree import Oncotree
from boostdm.globals import DISCOVERY_TIERS, MUTATION_TIERS, FSCORE_THRESHOLD


def get_fscore(model_evaluation):
    return np.nanmean(model_evaluation.get('fscore50', 0.5))


class Hierarchy:

    def __init__(self):
        self._tree = Oncotree()

    def climb(self, ttype, gene):
        """
        Args:
            ttype, gene
        Returns:
            generator of ascending {tumor-type, gene} pairs throughout the oncotree
            hierarchy starting from {ttype, gene}
        """

        while ttype is not None:
            yield ttype, gene
            ttype = self._tree.fetch_parent_ttype(ttype)


def meet_condition(fscore, discovery, n_muts):
    if fscore >= FSCORE_THRESHOLD:
        for discovery_thresh, n_muts_thresh in zip(DISCOVERY_TIERS, MUTATION_TIERS):
            if (discovery >= discovery_thresh) and (n_muts >= n_muts_thresh):
                return True
    return False


def evaluate(model_evaluations):

    res = {}

    hierarchy = Hierarchy()

    for ttype, gene in model_evaluations.keys():
        for tt, gg in hierarchy.climb(ttype, gene):

            fscore = model_evaluations[(tt, gg)]['fscore50']
            discovery = model_evaluations[(tt, gg)]['discovery']
            n_muts = model_evaluations[(tt, gg)]['n_muts']

            if meet_condition(fscore, discovery, n_muts):
                res[(ttype, gene)] = tt, gg
                break
    return res


@click.command()
@click.option('--eval_folder', 'eval_folder', type=click.Path(), help='input folder containing autoevaluation results')
@click.option('--discovery_path', 'discovery_path', type=click.Path(), help='file path to discovery output table')
@click.option('--output', 'output_file', type=click.Path(), help='output folder')
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

        models[(ttype, gene)] = {'fscore50': fscore50,
                                 'discovery': discovery,
                                 'n_muts': n_muts}

    res = evaluate(models)

    with gzip.open(output_file, 'wb') as f:
        pickle.dump(res, f)


if __name__ == '__main__':
    cli()
