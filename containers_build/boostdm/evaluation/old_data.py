import glob
import gzip
import os
import pickle

import click
import numpy as np
import pandas as pd

from boostdm.oncotree import Oncotree
from boostdm.globals import DRIVERS_PATH, DISCOVERY_TIERS, MUTATION_TIERS


def get_fscore(model_evaluation):
    return np.nanmean(model_evaluation.get('fscore50', 0.5))


class Hierarchy:

    def __init__(self):
        self._tree = Oncotree()
        drivers = pd.read_csv(DRIVERS_PATH, sep='\t')
        # self._roles = dict(zip(drivers.SYMBOL.values, drivers.ROLE.values))

    def climb(self, ttype, gene):
        """
        Args:
            ttype, gene
        Returns:
            generator of ascending {tumor-type, gene} pairs starting from {ttype, gene}
        """

        yield ttype, gene
        while ttype != 'CANCER':
            ttype = self._tree.fetch_parent_ttype(ttype)
            yield ttype, gene


def evaluate(model_evaluations, **params):

    res = {}

    hierarchy = Hierarchy()

    fscore_threshold = params['fscore_threshold']
    mutation_tiers = params['mutation_tiers']
    discovery_tiers = params['discovery_tiers']

    for ttype, gene in model_evaluations.keys():
        for tt, gg in hierarchy.climb(ttype, gene):
            fscore = model_evaluations[(tt, gg)]['fscore50']
            discovery = model_evaluations[(tt, gg)]['discovery']
            muts = model_evaluations[(tt, gg)]['n_muts']

            if fscore >= fscore_threshold:
                for disc_thresh, muts_thresh in zip(discovery_tiers, mutation_tiers):
                    if (discovery >= disc_thresh) and (muts >= muts_thresh):
                        res[(ttype, gene)] = tt, gg
                        break
        res[(ttype, gene)] = res.get((ttype, gene), (None, None))
    return res


@click.command()
@click.option('--eval_folder', 'eval_folder', type=click.Path(), help='input folder containing autoevaluation results')
@click.option('--discovery_path', 'discovery_path', type=click.Path(), help='file path to discovery output table')
@click.option('--fscore', type=float, help='F-score50 threshold')
@click.option('--output', 'output_file', type=click.Path(), help='output folder')
def cli(eval_folder, discovery_path, fscore, output_file):
    
    # discovery_tiers = (0, 0.2, 0.4, 0.6, 0.8)
    # mutation_tiers = (40, 30, 20, 10, 5)

    models = {}

    df_discovery = pd.read_csv(discovery_path, sep='\t')
    df_discovery = df_discovery.set_index(['gene', 'ttype'])

    for f in glob.glob(os.path.join(eval_folder, '*/*.eval.pickle.gz')):
        folder, file = os.path.split(f.replace(eval_folder, ''))
        ttype = folder.replace(os.path.sep, '')  # remove leading "/" if any
        gene = file.split('.')[0]
        with gzip.open(f, 'rb') as fd:
            d = pickle.load(fd)

        models[(ttype, gene)] = {'fscore50': get_fscore(d),
                                 'discovery': df_discovery.loc[(gene, ttype), 'discovery_index'],
                                 'n_muts': df_discovery.loc[(gene, ttype), 'n_muts']}

    res = evaluate(models,
                   fscore_threshold=fscore,
                   discovery_tiers=DISCOVERY_TIERS,
                   mutation_tiers=MUTATION_TIERS)

    with gzip.open(output_file, 'wb') as f:
        pickle.dump(res, f)


if __name__ == '__main__':
    cli()
