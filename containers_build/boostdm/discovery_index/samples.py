
import json
import operator
import functools
import click
import pandas as pd

from boostdm.oncotree import Oncotree


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

            if ttype not in samples_info:
                samples_info[ttype] = []

            list_hypermutators = []
            if 'hypermutators' in d[cohort]:
                list_hypermutators += list(d[cohort]['hypermutators']['hypermutators'])

            # samples_info[ttype] += list(d[cohort]['donors'].keys()) + list_hypermutators
            samples_info[ttype] += functools.reduce(operator.concat, d[cohort]['donors'].values()) + \
                                   list_hypermutators

    return samples_info


def aggregate_samples(samples):

    completion = samples
    tree = Oncotree()
    for ttype in tree.ttypes:
        if ttype not in samples:
            completion[ttype] = []
            for tt in tree.get_ttypes(ttype):
                completion[ttype] += samples.get(tt, [])
    return completion


@click.command()
@click.argument('files', nargs=-1)
@click.option('--cohorts', type=click.Path(exists=True))
@click.option('--output', type=click.Path())
def cli(cohorts, output, files):
    samples = build_table(cohorts, files)
    samples = aggregate_samples(samples)  # include more general tumor types
    with open(output, 'w') as fd:
        json.dump(samples, fd)


if __name__ == '__main__':
    cli()
