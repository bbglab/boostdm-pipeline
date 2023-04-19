
import json

import click
import pandas as pd


def build_table(cohorts, files):
    df_stats = pd.read_csv(cohorts, sep="\t")
    df_stats.rename(columns={"CANCER_TYPE": "cancer_type"}, inplace=True)

    # samples info
    samples = []
    for file in files:
        with open(file, 'r') as fd:
            samples.append(json.load(fd))

    samples_info = {}
    for d in samples:
        for cohort in d:
            try:
                ttype = df_stats[df_stats['COHORT'] == cohort]['cancer_type'].values[0]
            except:
                continue

            if not ttype in samples_info:
                samples_info[ttype] = []

            list_hypermutators = []
            if 'hypermutators' in d[cohort]:
                list_hypermutators += list(d[cohort]['hypermutators']['hypermutators'])
            samples_info[ttype] += list(d[cohort]['samples']['mut_per_sample'].keys()) + list_hypermutators
    return samples_info


@click.command()
@click.argument('files', nargs=-1)
@click.option('--cohorts', type=click.Path(exists=True))
@click.option('--output', type=click.Path())
def cli(cohorts, output, files):
    samples = build_table(cohorts, files)
    with open(output, 'w') as fd:
        json.dump(samples, fd)


if __name__ == '__main__':
    cli()
