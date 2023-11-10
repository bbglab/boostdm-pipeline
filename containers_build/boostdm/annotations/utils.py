
from boostdm.regression_utils import encode_feature


def encoding(df):
    """encoding categorical features for gradient boosting to handle them"""
    df = encode_feature(df, 'csqn_type')
    df = encode_feature(df, 'CLUSTL_cat')
    df = encode_feature(df, 'HotMaps_cat')
    df = encode_feature(df, 'smRegions_cat')
    # df = encode_feature(df, 'role')
    return df


def rectify_synonymous(df):

    dg = df.copy()
    ind = dg[dg['csqn_type_synonymous'] == 1].index
    forbidden = ['CLUSTL_SCORE', 'CLUSTL_cat_1',
                 'HotMaps_cat_1',
                 'smRegions_cat_1',
                 'nmd',
                 'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination']
                 # "cat_2" features for CLUSTL, Hotmaps and smRegions have been removed

    for c in forbidden:
        dg.loc[ind, c] = 0
    return dg


def rectify_missense(df):

    dg = df.copy()
    ind = dg[dg['csqn_type_missense'] == 1].index
    forbidden = ['nmd']
    for c in forbidden:
        dg.loc[ind, c] = 0
    return dg


def rectify_splicing(df):

    dg = df.copy()
    ind = dg[dg['csqn_type_splicing'] == 1].index
    forbidden = ['HotMaps_cat_1',
                 'nmd',
                 'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination']
                # "cat_2" features have been removed for Hotmaps
    for c in forbidden:
        dg.loc[ind, c] = 0
    return dg
