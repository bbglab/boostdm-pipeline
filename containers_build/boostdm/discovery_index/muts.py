
from os import path

import click
import pandas as pd


def build_table(files):
    mutations_table = []
    for file in files:
        cohort_name = path.basename(file).split('.')[0]
        df = pd.read_csv(file, sep='\t')
        df['COHORT'] = cohort_name
        mutations_table.append(df)
    df = pd.concat(mutations_table, axis=0)
    df['chr'] = df['chr'].astype(str)
    df['pos'] = df['pos'].astype(int)
    return df


@click.command()
@click.argument('files', nargs=-1)
@click.option('--output', type=click.Path())
def cli(files, output):
    df = build_table(files)
    df.to_csv(output, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    cli()
