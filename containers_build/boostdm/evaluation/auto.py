import gzip
import pickle

import click
import numpy as np
from sklearn.metrics import matthews_corrcoef, roc_auc_score, log_loss


def safe_roc_auc_score(x, y):
    """Safe version of roc_auc_score"""
    try:
        return roc_auc_score(x, y)
    except ValueError:
        return np.nan


def safe_log_loss(x, y):
    """Safe version of roc_auc_score"""
    try:
        return log_loss(x, y)
    except ValueError:
        return np.nan


def safe_accuracy(y_true, y_pred):

    tp = np.nansum(y_true * y_pred)
    fp = np.nansum((1 - y_true) * y_pred)
    tn = np.nansum((1 - y_true) * (1 - y_pred))
    fn = np.nansum(y_true * (1 - y_pred))
    try:
        return (tp + tn) / (tp + tn + fp + fn)
    except ValueError:
        return np.nan


def safe_precision(y_true, y_pred):

    tp = np.nansum(y_true * y_pred)
    fp = np.nansum((1 - y_true) * y_pred)
    try:
        return tp / (tp + fp)
    except ValueError:
        return np.nan


def safe_npv(y_true, y_pred):

    tn = np.nansum((1 - y_true) * (1 - y_pred))
    fn = np.nansum(y_true * (1 - y_pred))
    try:
        return tn / (tn + fn)
    except ValueError:
        return np.nan


def safe_recall(y_true, y_pred):

    tp = np.nansum(y_true * y_pred)
    fn = np.nansum(y_true * (1 - y_pred))
    try:
        return tp / (tp + fn)
    except ValueError:
        return np.nan


def safe_fscore(y_true, y_pred, beta=1.):

    precision = safe_precision(y_true, y_pred)
    recall = safe_recall(y_true, y_pred)
    try:
        return (1 + beta ** 2) * precision * recall * (1 / (precision * beta ** 2 + recall))
    except ValueError:
        return np.nan


def evaluate(model):
    """Compute several metrics of the models performances"""

    classifiers = model['models']
    covariates = model['x_test']
    labels = model['y_test']
    predicts = [model.predict_proba(x)[:, 1] for model, x in zip(classifiers, covariates)]

    res = dict()

    # ROC-AUC score
    res['auc'] = [safe_roc_auc_score(lab, pred) for lab, pred in zip(labels, predicts)]
    # MCC score
    res['mcc'] = [matthews_corrcoef(lab >= 0.5, pred >= 0.5) for lab, pred in zip(labels, predicts)]
    # Log-loss
    res['logloss'] = [safe_log_loss(lab, pred) for lab, pred in zip(labels, predicts)]
    # precision
    res['precision'] = [safe_precision(lab >= 0.5, pred >= 0.5) for lab, pred in zip(labels, predicts)]
    # NPV
    res['npv'] = [safe_npv(lab >= 0.5, pred >= 0.5) for lab, pred in zip(labels, predicts)]
    # recall
    res['recall'] = [safe_recall(lab >= 0.5, pred >= 0.5) for lab, pred in zip(labels, predicts)]
    # fscore
    res['fscore100'] = [safe_fscore(lab >= 0.5, pred >= 0.5) for lab, pred in zip(labels, predicts)]
    # fscore: beta=0.5
    res['fscore50'] = [safe_fscore(lab >= 0.5, pred >= 0.5, beta=0.5) for lab, pred in zip(labels, predicts)]
    # accuracy
    res['accuracy'] = [safe_accuracy(lab >= 0.5, pred >= 0.5) for lab, pred in zip(labels, predicts)]
    # test dataset balance
    res['balance'] = [0.5 - lab.mean() for lab, pred in zip(labels, predicts)]
    # calibration
    res['calibration'] = [(pred.mean() - lab.mean()) / lab.mean() for lab, pred in zip(labels, predicts)]
    # test size
    res['size'] = [len(y) for y in labels]

    return res


@click.command()
@click.option('--model', 'model_path', type=click.Path(exists=True), help='File corresponding to models for oncotree instance')
@click.option('--output', 'output_path', type=click.Path(), help='output file')
def cli(model_path, output_path):

    with gzip.open(model_path, 'rb') as f:
        model = pickle.load(f)

    res_dict = evaluate(model)

    with gzip.open(output_path, 'wb') as f:
        pickle.dump(res_dict, f)


if __name__ == '__main__':
    cli()
