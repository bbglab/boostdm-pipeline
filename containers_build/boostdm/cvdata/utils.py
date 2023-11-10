
import pandas as pd
import numpy as np

from boostdm.globals import COLUMNS_TRAINING


def sort_filter(x_train, x_test, y_train, y_test):
    """
    Final preparation steps to achieve a competent CV data:
    1) remove repeated data items in that test datasets
       removing duplicate sites in test dataset --but not in training-- as repeated data in training provide us with
       weight of evidence, whereas too many repeated data at testing can spoil our capacity to generalize well
    2) random sampling to get a balanced test set
    3) set training feature labels in a canonical order taken from configuration
    """

    # reset index

    x_test.reset_index(inplace=True, drop=True)
    y_test.reset_index(inplace=True, drop=True)

    # remove duplicates from test set

    x_test['chr'] = x_test['chr'].astype(str)
    x_test['pos'] = x_test['pos'].astype(int)
    x_test = x_test.drop_duplicates(['pos', 'alt'])
    test_index = x_test.index
    y_test = y_test.loc[y_test.index.intersection(test_index)]

    # balance test set

    total_index = set(x_test.index)
    if y_test.mean() <= 0.5:
        balance_index = y_test[y_test == 1].index.tolist()
    else:
        balance_index = y_test[y_test == 0].index.tolist()
    remaining_index = list(set(total_index) - set(balance_index))
    balance_index += list(np.random.choice(remaining_index, size=len(balance_index), replace=False))
    x_test = x_test.loc[x_test.index.intersection(balance_index)]
    y_test = y_test.loc[y_test.index.intersection(balance_index)]

    # feature labels in standard order
    avoid = ['chr', 'pos', 'ref', 'alt']
    features = list(filter(lambda x: x not in avoid, x_train))
    # assert (set(features) == set(COLUMNS_TRAINING))
    x_train = x_train[avoid + COLUMNS_TRAINING]
    x_test = x_test[avoid + COLUMNS_TRAINING]

    return x_train, x_test, y_train, y_test


def tuple_join(t, s):

    out = []
    assert (len(t) == len(s))
    for left, right in zip(t, s):
        df = pd.concat([left, right], sort=True, axis=0)
        df = df.reset_index(drop=True)
        out.append(df)
    return tuple(out)


def vertical_join(cvlist1, cvlist2):
    outcvlist = []
    for left, right in zip(cvlist1, cvlist2):
        outcvlist.append(tuple_join(left, right))
    return outcvlist
