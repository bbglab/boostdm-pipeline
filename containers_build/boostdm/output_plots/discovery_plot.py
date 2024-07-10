import os
import json
import click

import glob
from functools import reduce
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

from boostdm.output_plots.utils.discovery_utils import plot_fit

# graphical config

def config_params(font_size=12):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

# cli

@click.command()
@click.option('--gene', type=str)
@click.option('--ttmodel', type=str)
def cli(gene, ttmodel):

    config_params()

    # samples

    with open(os.path.join(os.environ['OUTPUT'], 'discovery', 'samples.json'), 'rt') as g:
        samples = json.load(g)

    # mutations

    mutations = pd.read_csv(os.path.join(os.environ['OUTPUT'], 'discovery', 'mutations.tsv.gz'), sep='\t')

    # plot

    fig, ax = plt.subplots(figsize=(3, 3))
    plot_fit(gene, ttmodel, samples, mutations, ax)
    plt.savefig(f'{gene}.{ttmodel}.bending.svg', bbox_inches='tight', dpi=300)
    plt.savefig(f'{gene}.{ttmodel}.bending.png', bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == '__main__':
    
    cli()
