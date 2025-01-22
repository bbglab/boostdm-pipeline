from collections import defaultdict
from os import path

import click
import pandas as pd

from boostdm.features import clustl, hotmaps, smregions

FEATURES_DICT = {
    'clustl' : clustl,
    'hotmaps': hotmaps,
    'smregions' : smregions
}


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

@click.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.05, type=float, help='Pvalue/Qvalue threshold')
@click.option('--cohorts', type=click.Path(), required=True, help='cohorts file')
@click.option(
    '--method',
    type=click.Choice(list(FEATURES_DICT.keys()), case_sensitive=False),
    required=True,
    help=f'Methods used to compute features. Choose one of: {FEATURES_DICT.keys()}'
)
@click.argument('files', nargs=-1)
def cli(files, output, threshold, cohorts, method):
    """Group features for either ClustL, Hotmaps and SMRegions."""

    feature = FEATURES_DICT.get(method)

    data = []
    ttypes_map = load_ttypes_map(cohorts)
    for ttype, files_list in group_by_ttype(files, ttypes_map):

        df = feature.generate(files_list, thresh=threshold )
        df['CANCER_TYPE'] = ttype
        data.append(df)

    df = pd.concat(data, axis=0)
    df.to_csv(output, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
