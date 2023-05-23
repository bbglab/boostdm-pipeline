
import click
import pandas as pd


def matching_cluster(clusts, chr_, pos):
    """
    clusts: dataframe
    chr_: chromosome
    pos: position
    returns a tuple (bool, value):
        bool: whether there is a cluster overlapping the genomic coords (chr_, pos)
        value: if the bool is True, value = clusters overlapping the genomic coords
               if the bool is False, value = None
    """

    # check for chr_ and position in cluster table
    if (clusts is not None) and (len(clusts) > 0):
        c = clusts[(clusts['CHROMOSOME'] == chr_) & (clusts['5_COORD'] <= pos) & (pos <= clusts['3_COORD'])]
        if len(c) > 0:
            return True, c
    return False, None


def assign(chr_, pos, clusters):

    e1, _ = matching_cluster(clusters['specific'], chr_, pos)
    e2, _ = matching_cluster(clusters['pan'], chr_, pos)


    # debugging

    if pos == 7829993:

        print(matching_cluster(clusters['specific'], chr_, pos))
        print(matching_cluster(clusters['pan'], chr_, pos))

    if e1:
        return 2.  # cluster in specific tumor type
    elif e2:
        return 1.  # cluster in another tumor type
    else:
        return 0.  # no cluster


def add_feature(df, specific_df, global_df):

    clusters = {
        'specific': specific_df,
        'pan': global_df
    }

    df['CLUSTL'] = df.apply(
        lambda x: pd.Series(assign(x['chr'], x['pos'], clusters)), 
        axis=1)
    
    return df


def generate(files, pval_thresh=0.05):
    """group a set of output files into a single dataframe"""

    df = pd.DataFrame()

    for file in files:
        input1 = pd.read_csv(file,
                             sep='\t',
                             usecols=['CHROMOSOME', 'COORDINATES', 'SCORE', 'P']).drop_duplicates()
        if input1.shape[0] == 0:
            continue
        input2 = input1.COORDINATES.str.split(';').apply(pd.Series)
        input2.index = input1.set_index(['CHROMOSOME', 'SCORE', 'P']).index
        input3 = input2.stack().reset_index(['CHROMOSOME', 'SCORE', 'P'])
        input3['5_COORD'], input3['3_COORD'] = zip(*input3[0].map(lambda x: tuple(x.split(','))))
        del input3[0]
        df = df.append(input3[input3.P <= pval_thresh], ignore_index=True, sort=False)

    if len(df) == 0:
        return pd.DataFrame(columns=['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE'])
    else:
        df = df.sort_values(by=['CHROMOSOME', '5_COORD'])
        df['5_COORD'] = df['5_COORD'].astype('int32')
        df['3_COORD'] = df['3_COORD'].astype('int32')
        df['CHROMOSOME'] = df['CHROMOSOME'].apply(str)
        return df[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]


@click.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.05, type=float, help='Pvalue threshold')
@click.argument('files', nargs=-1)
def cli(files, output, threshold):

    df = generate(files, pval_thresh=threshold)

    df.to_csv(output, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
