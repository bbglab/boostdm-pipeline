
import json
import os
from functools import partial

import numpy as np
import pandas as pd

from boostdm.globals import PTMS_FILE

with open(PTMS_FILE, 'rt') as read_file:
    DATA = json.load(read_file)


def get_ptms(d_data, hugo, mutation, restrict_missense=False):

    """
    Given a HUGO and a mutation returns the Post Translational Modification Sistes (PTMs) affected by the current mutation.
    Specific alterations considered: 'Acetylation', 'O-glycosylation', 'O-GlcNAc', 'Phosphorylation',
                                     'Ubiquitination', 'Methylation', 'Sumoylation'
    One global column called 'Regulatory_Site' encompasses other type of regulatory sites.

    Params:
    -------
    hugo: HUGO symbol of the query gene
    mutation: Mutation with the format A33C. The WT and MT amino acids are required. Synonymous mutations are encoded as A33A.
              Nosense mutations and frameshift are encoded as A33* or A33-
    restrict_missense: whether the search should be restricted to missense mutations or all type of mutations

    Returns:
    --------
    List of affected PTMs in the following order ['Acetylation', 'Phosphorylation', 'Ubiquitination',
                                                  'Methylation', 'Regulatory_Site'].
    Values: affected=1; non-affected=0
    """

    if mutation is None:
        return np.zeros(5)

    aa_wt = mutation[0]
    aa_mt = mutation[-1]
    aa = mutation[:-1]
    if restrict_missense:
        if (aa_mt == aa_wt) or (aa_mt == "*") or (aa_mt == "-"):
            return np.zeros(5)
    if (hugo in d_data) and (aa in d_data[hugo]):
        data_vector = np.array(d_data[hugo][aa])
        return data_vector[[0, 5, 3, 7, 4]]  # 0: acetylation; 5: methylation; 3: phosphorylation;
                                             # 7: regulatory site; # 4: ubiquitination;
    else:
        return np.zeros(5)


def add_ptms(df):

    ptms_cols = ['Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination']

    get_annotation = partial(get_ptms, DATA)
    df[ptms_cols] = df[['gene', 'aachange']].apply(lambda r: pd.Series(get_annotation(*r)), axis=1)
    return df


if __name__ == '__main__':

    df = pd.DataFrame({'gene': ['ABCA10'], 'aachange': ['G539W']})
    df = add_ptms(df)
    print(df)

