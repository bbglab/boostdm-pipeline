import pandas as pd
import numpy as np

from src.globals import COLUMNS_TRAINING


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
    x_test = x_test.reset_index(drop=True)
    y_test = y_test.reset_index(drop=True)

    # remove duplicates from test set

    x_test["chr"] = x_test["chr"].astype(str)
    x_test["pos"] = x_test["pos"].astype(int)
    x_test = x_test.drop_duplicates(subset=["pos", "alt"], keep="first")
    y_test = y_test.loc[x_test.index]

    # Balance the test set
    positive_idx = y_test[y_test == 1].index
    negative_idx = y_test[y_test == 0].index

    # Find the minority class size
    min_size = min(len(positive_idx), len(negative_idx))

    # Downsample the majority class to match the minority size
    if len(positive_idx) > min_size:
        positive_idx = np.random.choice(positive_idx, size=min_size, replace=False)
    else:
        negative_idx = np.random.choice(negative_idx, size=min_size, replace=False)

    # Combine balanced indices
    balanced_idx = np.concatenate([positive_idx, negative_idx])

    # Apply the balanced index selection
    x_test = x_test.loc[balanced_idx].copy()
    y_test = y_test.loc[balanced_idx].copy()

    # feature labels in standard order
    avoid = ["chr", "pos", "ref", "alt"]
    x_train = x_train[avoid + COLUMNS_TRAINING].copy()
    x_test = x_test[avoid + COLUMNS_TRAINING].copy()

    return x_train, x_test, y_train, y_test


def tuple_join(t, s):
    out = []
    assert len(t) == len(s)
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
