import os
import glob
import click
import json

import numpy as np
import pandas as pd

from boostdm.globals import INTOGEN_DATASETS, DRIVERS_PATH


drivers = pd.read_csv(DRIVERS_PATH, sep='\t')
driver_genes = drivers['SYMBOL'].unique()

csqn_type_dict = {
    'Missense': 'missense',
    'Synonymous': 'synonymous',
    'Nonsense': 'nonsense',
    'Essential_Splice': 'splicing',
    'Stop_loss': 'nonsense',
    'no-SNV': 'non_snv'
}

@click.command()
@click.option('--output', type=click.Path())
@click.option('--percentile', type=int)
def cli(output, percentile=5):

    # retrieve mutations

    total_df = []
    for fn in glob.glob(f"{INTOGEN_DATASETS}/steps/dndscv/*.dndscv_annotmuts.tsv.gz"):
        df = pd.read_csv(fn, sep='\t')
        total_df.append(df)
    total_df = pd.concat(total_df, axis=0)
    total_df['impact'] = total_df['impact'].apply(lambda x: csqn_type_dict[x])
    total_df['pos'] = total_df['pos'].astype(int)
    total_df['chr'] = total_df['chr'].astype(str)
    total_df = total_df[total_df['gene'].isin(driver_genes)]  # intogen driver genes only
    total_snvs = total_df[total_df['impact'] != 'non_snv']  # snvs only

    # mutation count dict

    d = total_snvs.groupby(['gene', 'impact']).size().to_dict()

    gene_dict = {}  # dict: gene -> dict csqn_type -> mut count
    for k, v in d.items():
        g = k[0]
        gene_dict[g] = {**gene_dict.get(g, {}), **{k[1]: v}}

    # compute delta log counts and TMBs

    y, z = [], []
    tmby, tmbz = [], []
    for g, v in gene_dict.items():
        deltanon = np.log10(v.get('nonsense', 1)) - np.log10(v['synonymous'])
        deltaspl = np.log10(v.get('splicing', 1)) - np.log10(v['synonymous'])
        y.append(deltanon), z.append(deltaspl)

    thresh_non, thresh_spl = np.quantile(y, percentile / 100), np.quantile(z, percentile / 100)

    csqn_type_vetting = {}
    for g, v in gene_dict.items():
        csqn_type_vetting[g] = []
        deltanon = np.log10(v.get('nonsense', 1)) - np.log10(v['synonymous'])
        lognon = np.log10(v.get('nonsense', 1))
        if (deltanon <= thresh_non) or (lognon == 0):
            csqn_type_vetting[g].append('nonsense')
        deltaspl = np.log10(v.get('splicing', 1)) - np.log10(v['synonymous'])
        logspl = np.log10(v.get('splicing', 1))
        if (deltaspl <= thresh_spl) or (logspl == 0):
            csqn_type_vetting[g].append('splicing')

    with open(output, 'w') as f:
        json.dump(csqn_type_vetting, f)

if __name__ == '__main__':
    
    cli()
