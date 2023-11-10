"""Script to conduct training of gradient boosting models"""


# Imports

import gzip
import os
from multiprocessing import Pool

import click
import dill as pickle
import numpy as np
import xgboost as xgb

import warnings
warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
warnings.filterwarnings(module='sklearn*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='pandas*', action='ignore', category=RuntimeWarning)

from boostdm.globals import XGB_PARAMS

# TODO: get rid of this dependency using xgboost directly!
from boostwrap import Classifier


# Globals
# -------

non_features = ['pos', 'chr', 'ref', 'alt']

# Classification Training Functions
# ---------------------------------


def train(values):
    """Returns: optimal classification threshold and trained XGB model"""

    x_train, x_test, y_train, y_test, split_number, seed = tuple(values)

    if x_test.shape[0] == 0:
        return None

    params = XGB_PARAMS.copy()
    params['n_estimators'] = 20000  # set it high enough to allow "early stopping" events below
    params['base_score'] = y_train.mean()
    params['n_jobs'] = 1
    params['seed'] = seed
    myclassifier = Classifier(**params)

    # train with xgboost
    learning_curve_dict = {}
    myclassifier.train(x_train, y_train,
                       eval_set=[(x_train, y_train), (x_test, y_test)],
                       eval_metric='logloss',  # mcc_loss could be used here
                       early_stopping_rounds=2000,
                       callbacks=[
                           xgb.callback.record_evaluation(learning_curve_dict),
                       ],
                       verbose=False)

    params['n_estimators'] = myclassifier.model.best_iteration
    learning_curve_dict = {k: v['logloss'][:params['n_estimators']] for k, v in learning_curve_dict.items()}
    myclassifier.model.set_params(**params)

    return myclassifier, split_number, x_test, y_test, learning_curve_dict


# Interactive Commands: training, bgqmap_training
# -----------------------------------------------

@click.command()
@click.option('--splits', 'file_cv', type=click.Path(exists=True), help='Folder with cvdata', required=True)
@click.option('--output', 'output_file', type=click.Path(), help='File name for the models', required=True)
@click.option('--cores', type=int, help='Number of cores to be used', default=1)
@click.option('--min-rows', type=int, help='Minimum number of rows to carry out training', default=30)
@click.option('--seed', type=int, default=None)
def cli(file_cv, output_file, cores, min_rows, seed):

    np.random.seed(seed)

    if not os.path.exists(file_cv):
        raise FileExistsError('Missing splits file')

    dict_results = {'models': [], 'split_number': [], 'x_test': [], 'y_test': [], 'learning_curves': []}

    # List the genes used in the computation
    with Pool(cores) as p:

        # Load the file after the pool is created so it gets not replicated in the processes
        with gzip.open(file_cv, 'rb') as f:
            split_cv = pickle.load(f)

        if split_cv is None:
            raise Exception('Cross-validation splits are not defined (None)')

        mean_size = np.nanmean([cv[0].shape[0] for cv in split_cv])
        print(mean_size, min_rows)

        if mean_size < min_rows:
            raise Exception('Mean training size is not sufficiently high')

        list_cvs = []
        for i, x in enumerate(split_cv):
            x_list = list(x) + [i, np.random.randint(100000)]

            # filter out non-features, i.e., columns not used for training
            x_list[0] = x_list[0].drop(non_features, axis=1)
            x_list[1] = x_list[1].drop(non_features, axis=1)
            list_cvs.append(x_list)

        for out in p.imap(train, list_cvs):
            if out is None:
                pass
            else:
                model, split_number, x_test, y_test, learning_curve = out
                dict_results['models'].append(model)
                dict_results['split_number'].append(split_number)
                dict_results['x_test'].append(x_test)
                dict_results['y_test'].append(y_test)
                dict_results['learning_curves'].append(learning_curve)

    with gzip.open(output_file, 'wb') as f:
        pickle.dump(dict_results, f)

if __name__ == '__main__':
    cli()
