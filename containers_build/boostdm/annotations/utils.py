
import pandas as pd


csqn_type_list = ['missense', 'nonsense', 'splicing', 'synonymous']


def encode_consequence_type(data):

    data.loc[~data['csqn_type'].isin(csqn_type_list), 'csqn_type'] = 'none'
    one_hot = pd.get_dummies(data, columns=['csqn_type'], prefix_sep='_')
    one_hot.drop(columns=['csqn_type_none'], inplace=True, errors='ignore')
    for c in csqn_type_list:
        col = f'csqn_type_{c}'
        if col not in list(one_hot.columns):
            one_hot[col] = 0
    return one_hot


def rectify_synonymous(df):

    ind = df.index[df['csqn_type_synonymous'] == 1].tolist()
    forbidden = [
        'CLUSTL', 'HotMaps', 'smRegions',
        'nmd',
        'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination'
        ]
    for c in forbidden:
        df.loc[ind, c] = 0
    return df


def rectify_missense(df):

    ind = df.index[df['csqn_type_missense'] == 1].tolist()
    forbidden = ['nmd']
    for c in forbidden:
        df.loc[ind, c] = 0
    return df


def rectify_splicing(df):

    ind = df.index[df['csqn_type_splicing'] == 1].tolist()
    forbidden = ['HotMaps', 'nmd', 'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination']
    for c in forbidden:
        df.loc[ind, c] = 0
    return df
