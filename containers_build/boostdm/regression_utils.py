import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.metrics import matthews_corrcoef, roc_curve, roc_auc_score
import warnings
warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
warnings.filterwarnings(module='sklearn*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='matplotlib*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='pandas*', action='ignore', category=RuntimeWarning)
warnings.filterwarnings(module='ipykernel*', action='ignore', category=FutureWarning)
warnings.filterwarnings(module='shap*', action='ignore', category=RuntimeWarning)


csqn_type_list = ['missense', 'nonsense', 'splicing', 'synonymous']


def encode_consequence_type(data, feature):
    
    data.loc[~data[feature].isin(csqn_type_list), feature] = 'none'
    one_hot = pd.get_dummies(data[[feature]], columns=[feature], prefix_sep='_')
    for value in csqn_type_list:
        col = feature + '_' + str(value)
        if col not in list(one_hot.columns):
            one_hot[col] = 0
    data.drop(columns=feature, inplace=True)
    to_remove = feature + '_' + 'none'
    if to_remove in one_hot.columns:
        one_hot.drop(columns=to_remove, inplace=True)
    canonical_order = [feature + '_' + str(value) for value in csqn_type_list]
    data = data.join(one_hot[canonical_order])
    return data


def encoding_test():

    data_dict = {'feat1': [1., 2., 3.], 'feat2': [1.4, 3.1, 3.2], 'csqn_type': ['strange', 'none', 'splice_acceptor_variant']}
    data = pd.DataFrame(data_dict)
    df = encode_consequence_type(data, 'csqn_type')
    print(df.columns)


def propagate_slice(x_data, y_data, sliced_data):
    """Given a subset of data, it slices x_data and y_data accordingly"""

    subset_index = [i for i in x_data.index if i in sliced_data.index]
    return x_data.loc[subset_index, :], y_data.loc[subset_index]


# Evaluation
# ----------


def matthews(tp, fp, tn, fn):
    return (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))


def accuracy(tp, fp, tn, fn):
    return (tp + tn) / (tp + fp + tn + fn)


# MCC Score

def mcc_score(model, x_test, y_test, thresh=0.5):
    dtest = xgb.DMatrix(x_test.values, label=y_test.values)
    preds = model.predict_proba(x_test)[:, 1]
    labels = dtest.get_label()
    return matthews_corrcoef(labels, preds > thresh)


# Evaluation: ROC Score

def roc_score(model, x_test, y_test):
    dtest = xgb.DMatrix(x_test.values, label=y_test.values)
    preds = model.predict_proba(x_test)[:, 1]
    labels = dtest.get_label()
    return roc_auc_score(labels, preds)


def roc_points(model, x_test, y_test):
    """Returns: fprs, tprs and thresholds"""

    dtest = xgb.DMatrix(x_test.values, label=y_test.values)
    preds = model.predict_proba(x_test)[:, 1]
    labels = dtest.get_label()
    return roc_curve(labels, preds)


def roc_plot(myclassifier, x_test, y_test, title=True):
    fprs, tprs, thresholds = roc_points(myclassifier, x_test, y_test)
    plt.plot(fprs, tprs)
    s = roc_score(myclassifier, x_test, y_test)
    if title:
        plt.title(f'ROC Curve: AUC = {np.round(s,2)}')
    plt.xlabel('FPR')
    plt.ylabel('TPR')


def compute_test_summary(myclassifier, x_test, y_test, thresh=0.5):
    cols = ['predicted_prob', 'predicted_binary', 'true_condition']
    test_result = pd.DataFrame(columns=cols)
    test_result['predicted_prob'] = myclassifier.predict_proba(x_test)[:, 1]
    test_result['predicted_binary_thresh'] = test_result['predicted_prob'].apply(lambda x: int(x >= thresh))
    test_result['predicted_binary_standard'] = myclassifier.predict(x_test)
    test_result['true_condition'] = y_test.values
    return test_result


def deviance_loss(test_summary):
    df = test_summary.copy()
    df['deviance'] = df.apply(lambda v: -2 * np.log(v['predicted_prob']) * (v['true_condition'] == 1) - 2 * np.log(1 - v['predicted_prob']) * (v['true_condition'] == 0), axis=1)
    return df['deviance'].mean()


def evaluation(myclassifier, x_test, y_test, thresh=0.5):
    """
    myclassifier (boostwrap.Classifier)
    x_test (pandas dataframe)
    y_test (pandas series)
    """

    test_summary = compute_test_summary(myclassifier, x_test, y_test, thresh=thresh)

    contingency = pd.crosstab(test_summary.true_condition, test_summary.predicted_binary_thresh)

    tp = contingency.loc[1, 1]
    fp = contingency.loc[0, 1]
    fn = contingency.loc[1, 0]
    tn = contingency.loc[0, 0]

    return (contingency,
            accuracy(tp, fp, tn, fn),
            matthews(tp, fp, tn, fn),
            deviance_loss(test_summary))


if __name__ == '__main__':
    encoding_test()
