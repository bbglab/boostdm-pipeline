from collections import defaultdict
from os import path

import click
import pandas as pd

from boostdm.features import clustl, hotmaps, smregions


@click.group()
def cli():
    pass


def load_ttypes_map(file):
    df = pd.read_csv(file, sep='\t')
    df = df[['COHORT', 'CANCER_TYPE']]
    return df.set_index('COHORT').to_dict()['CANCER_TYPE']


def group_by_ttype(files, ttype_map):
    groups = {}
    for file in files:
        cohort_name = path.basename(file).split('.')[0]
        ttype = ttype_map[cohort_name]
        groups[ttype] = groups.get(ttype, []) + [file]

    for group, files_list in groups.items():
        yield group, files_list


@cli.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.05, type=float, help='Pvalue threshold')
@click.option('--cohorts', type=click.Path(), required=True, help='cohorts file')
@click.argument('files', nargs=-1)
def group_clustl(files, output, threshold, cohorts):

    data = []
    ttypes_map = load_ttypes_map(cohorts)
    for ttype, files_list in group_by_ttype(files, ttypes_map):

        df = clustl.generate(files_list, pval_thresh=threshold)
        df['CANCER_TYPE'] = ttype
        data.append(df)

    df = pd.concat(data, axis=0)
    df.to_csv(output, sep='\t', index=False, compression="gzip")


@cli.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.05, type=float, help='Pvalue threshold')
@click.option('--cohorts', type=click.Path(), required=True, help='cohorts file')
@click.argument('files', nargs=-1)
def group_hotmaps(files, output, threshold, cohorts):

    data = []
    ttypes_map = load_ttypes_map(cohorts)
    for ttype, files_list in group_by_ttype(files, ttypes_map):
        df = hotmaps.generate(files_list, pval_thresh=threshold)
        df['CANCER_TYPE'] = ttype
        data.append(df)

    df = pd.concat(data, axis=0)
    df.to_csv(output, sep='\t', index=False, compression="gzip")


@cli.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.05, type=float, help='Pvalue threshold')
@click.option('--cohorts', type=click.Path(), required=True, help='cohorts file')
@click.argument('files', nargs=-1)
def group_smregions(files, output, threshold, cohorts):

    data = []
    ttypes_map = load_ttypes_map(cohorts)
    for ttype, files_list in group_by_ttype(files, ttypes_map):
        df = smregions.generate(files_list, qval_thresh=threshold)
        df['CANCER_TYPE'] = ttype
        data.append(df)

    df = pd.concat(data)
    df.to_csv(output, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
