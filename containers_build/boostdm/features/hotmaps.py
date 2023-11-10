import click
import pandas as pd


def assign(chr_, pos, clusters):

    clusters_specific = clusters['specific']
    cluster1 = clusters_specific[
        (clusters_specific['chromosome'] == chr_) & (clusters_specific['pos'] == pos)]

    if len(cluster1) > 0:
        return 1
    else:
        clusters_pan = clusters['pan']
        cluster2 = clusters_pan[
            (clusters_pan['chromosome'] == chr_) & (clusters_pan['pos'] == pos)]
        if len(cluster2) > 0:
            return 2

    return 3  # no cluster


def add_feature(df, specific_df, global_df):

    clusters = {
        'specific': specific_df,
        'pan': global_df
    }

    df['HotMaps_cat'] = df.apply(
        lambda x: pd.Series(assign(x['chr'], x['pos'], clusters)), axis=1
    )

    df.HotMaps_cat = df.HotMaps_cat.astype('int8')

    return df


def generate(files, pval_thresh=0.05):
    
    """group a set of output files into a single dataframe"""
    df = pd.DataFrame()

    for input_file in files:
        input = pd.read_csv(input_file, sep='\t', header=0,
                            usecols=['chromosome', 'genomic position', 'q-value']).drop_duplicates()
        df = df.append(input[input['q-value'] <= pval_thresh], ignore_index=True, sort=False)

    if df.empty:
        df = pd.DataFrame(columns=['chromosome', 'pos'])
    else:
        s1 = df['genomic position'].str.split(',', expand=True).stack().str.strip().reset_index(level=1, drop=True)
        df1 = pd.DataFrame(s1, columns=['pos'])
        df = df.drop(['genomic position'], axis=1).join(df1).reset_index(drop=True)
        df.chromosome = df.chromosome.str.replace('chr', '')
        df.pos = df.pos.astype('int32')
        df = df.sort_values(by=['chromosome', 'pos'])

    return df[['chromosome', 'pos']]


@click.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.05, type=float, help='Pvalue threshold')
@click.argument('files', nargs=-1)
def cli(files, output, threshold):

    df = generate(files, pval_thresh=threshold)

    df.to_csv(output, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
