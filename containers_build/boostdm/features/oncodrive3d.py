import click
import pandas as pd


def assign(chr_, pos, clusters):

    clusters_specific = clusters['specific']
    cluster1 = clusters_specific[(clusters_specific['chromosome'] == chr_) & (clusters_specific['pos'] == pos)]

    if len(cluster1) > 0:
        return 2.  # tumor type specific cluster
    
    else:
        clusters_pan = clusters['pan']
        cluster2 = clusters_pan[(clusters_pan['chromosome'] == chr_) & (clusters_pan['pos'] == pos)]
        if len(cluster2) > 0:
            return 1.  # another tumor type cluster

    return 0.  # no cluster


def add_feature(df, specific_df, global_df):

    clusters = {
        'specific': specific_df,
        'pan': global_df
    }

    df['Oncodrive3D'] = df.apply(
        lambda x: pd.Series(assign(x['chr'], x['pos'], clusters)), 
        axis=1)

    return df


def generate(files, pval_thresh=0.05):
    """group a set of output files into a single dataframe"""
    df = pd.DataFrame()

    for input_file in files:
        input = pd.read_csv(input_file, sep='\t', header=0, usecols=['Gene', 'C_gene', 'C_pos']).drop_duplicates()
        df = df.append(input[input['C_gene'] == 1], ignore_index=True, sort=False)

    if df.empty:
        df = pd.DataFrame(columns=['Gene', 'C_gene', 'C_pos'])

    return df[['Gene', 'C_gene', 'C_pos']]