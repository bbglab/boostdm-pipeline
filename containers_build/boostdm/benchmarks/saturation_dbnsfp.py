# generate saturation mutagenesis with dbNSFP annotations

import tqdm
import click
import os
import glob

import numpy as np
import pandas as pd

from boostdm import BoostDMError


bad_format_scores = ['AlphaMissense_score', 'ESM1b_score', 'EVE_score', 'FATHMM_score', 'MetaRNN_score', 
                     'MutationAssessor_score', 'PROVEAN_score', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 
                     'REVEL_score', 'SIFT4G_score', 'SIFT_score', 'VEST4_score']

basepath = os.environ['OUTPUT']


def mean_dbNSFP_transform(r):
    """
    Mean summarization for dbNSFP columns in VEP output
    """
    
    if r != '-':
        values = [float(c) for c in r.split(',') if c != '.']
        if len(values) > 0:
            return np.mean(values)
    return np.nan


@click.command()
@click.option('--input', type=click.Path())
def cli(input):

    gene = os.path.basename(input).split('.')[0]
    pool = []
    for fn in glob.glob(os.path.join(basepath, f'saturation/prediction/{gene}.model.*.features.*.prediction.tsv.gz')):
            ttype_model = os.path.basename(fn).split('.')[2]
            ttype_features = os.path.basename(fn).split('.')[4]
            pool.append((gene, ttype_model, ttype_features))

    if len(pool) == 0:
        print(f'The saturation pool for gene={gene} is empty')
        return
    
    for gene, ttype_model, _ in tqdm.tqdm(pool):

        saturation_annotation_fn = os.path.join(basepath, 'saturation', 'annotation', f'{gene}.{ttype_model}.annotated.tsv.gz')
        sat_df = pd.read_csv(saturation_annotation_fn, sep='\t')

        annot_df = pd.read_csv(input, sep='\t')

        # reformatting

        annot_df.reset_index(inplace=True, drop=True)
        annot_df['chr'], annot_df['pos'] = zip(*annot_df['Location'].apply(lambda r: r.split(':')))
        annot_df = annot_df.rename(columns={'Allele': 'alt'})
        annot_df['chr'] = annot_df['chr'].astype(str)
        annot_df['pos'] = annot_df['pos'].astype(int)

        sat_df.reset_index(inplace=True, drop=True)
        sat_df['chr'] = sat_df['chr'].astype(str)
        sat_df['pos'] = sat_df['pos'].astype(int)

        # filter MANE_SELECT

        annot_df = annot_df[annot_df['MANE_SELECT'] != '-']

        # merge with dbNSFP

        sat_annot_df = pd.merge(sat_df, annot_df, on=['chr', 'pos', 'alt'], suffixes=('', '_y'))
        sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)

        for c in bad_format_scores:
            sat_annot_df[c] = sat_annot_df[c].apply(mean_dbNSFP_transform)
        sat_annot_df = sat_annot_df.replace('-', np.nan)

        sat_annot_df.to_csv(f'{gene}.{ttype_model}.saturation.dbNSFP.tsv', sep='\t', index=False)


if __name__ == '__main__':

    cli()
