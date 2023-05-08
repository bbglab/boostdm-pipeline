
from boostdm.regression_utils import encode_cluster, encode_consequence_type


def encoding(df):
    """encoding categorical features for gradient boosting to handle them"""
    df = encode_consequence_type(df, 'csqn_type')
    df = encode_cluster(df, 'CLUSTL_cat')
    df = encode_cluster(df, 'HotMaps_cat')
    df = encode_cluster(df, 'smRegions_cat')
    return df


def rectify_synonymous(df):

    dg = df.copy()
    ind = dg[dg['csqn_type_synonymous'] == 1].index
    forbidden = ['CLUSTL_SCORE', 'CLUSTL_cat_1', 'CLUSTL_cat_2',
                 'HotMaps_cat_1', 'HotMaps_cat_2',
                 'smRegions_cat_1', 'smRegions_cat_2',
                 'nmd',
                 'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination']
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
    forbidden = ['HotMaps_cat_1', 'HotMaps_cat_2',
                 'nmd',
                 'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination']
    for c in forbidden:
        dg.loc[ind, c] = 0
    return dg
