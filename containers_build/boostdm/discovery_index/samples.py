
import json
import operator
import functools
import click
import pandas as pd

from boostdm.oncotree import Oncotree
from boostdm.globals import COHORTS_PATH


def target_ttypes(cohort):
    oncotree = Oncotree()
    res = []
    for ttype in oncotree.ttypes:
        ttype_cohorts = oncotree.get_cohorts(ttype)
        if cohort in ttype_cohorts:
            res.append(ttype)
    return res


def build_table(input_variants):

    cohorts = pd.read_csv(COHORTS_PATH, sep='\t')
    cohorts_ttype_dict = dict(zip(cohorts['COHORT'], cohorts['CANCER_TYPE']))

    # samples info
    with open(input_variants, 'r') as fd:
        samples = json.load(fd)
    
    samples_info = {}
    for cohort in samples:

        # the samples of the cohort have to be appended 
        # to all the tumor types where the cohorts belongs

        for ttype in target_ttypes(cohort):
            if ttype not in samples_info:
                samples_info[ttype] = []
            list_hypermutators = []
            if 'hypermutators' in samples[cohort]:
                list_hypermutators += list(samples[cohort]['hypermutators']['hypermutators'])
            samples_info[ttype] += functools.reduce(operator.concat, samples[cohort]['donors'].values()) + list_hypermutators
    
    return samples_info


@click.command()
@click.option('--input', type=click.Path(exists=True))
@click.option('--output', type=click.Path())
def cli(input, output):
    samples = build_table(input)
    with open(output, 'w') as fd:
        json.dump(samples, fd)


if __name__ == '__main__':
    cli()
