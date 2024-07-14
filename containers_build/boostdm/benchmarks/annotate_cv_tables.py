import os
import glob
import tqdm

import pandas as pd


tp53_assays = ['TP53_Boettcher', 'TP53_Giacomelli', 'TP53_Kato', 'TP53_Kotler', 'TP53_Ursu']
dnmt3a_assays = ['DNMT3A_Lue']
ras_assays = ['Ras_Bandaru', 'KRAS_Ursu']
pten_assays = ['PTEN_Mighell', 'PTEN_Matreyek']


all_dbnsfp_scores = ['AlphaMissense_score', 'CADD_raw', 'ESM1b_score', 'EVE_score', 'FATHMM_score', 'MetaLR_score',
                     'MetaRNN_score', 'MutationAssessor_score', 'PROVEAN_score', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 
                     'REVEL_score', 'SIFT4G_score', 'SIFT_score', 'VEST4_score', 
                     'phyloP100way_vertebrate', 'phyloP17way_primate', 'phyloP470way_mammalian']


def cli():

    for fn in tqdm.tqdm(glob.glob(os.path.join(os.environ['OUTPUT'], 'benchmarks', 'cv_tables', '*.*.50.iter.tsv'))):
        
        gene, ttype = os.path.basename(fn).split('.')[:2]

        # cross-validation data

        cv_df = pd.read_csv(fn, sep='\t', low_memory=False)
        cv_df['chr'] = cv_df['chr'].astype(str)
        
        # annotate MAVE data

        saturation_mave_fn = os.path.join(os.environ['OUTPUT'], 'benchmarks', 'saturation_mave', f'{gene}.{ttype}.saturation.mave.tsv')
        if os.path.isfile(saturation_mave_fn):
            
            sat_annot = pd.read_csv(saturation_mave_fn, sep='\t')
            sat_annot['chr'] = sat_annot['chr'].astype(str)
                        
            if gene == 'TP53':
                cv_annot = pd.merge(cv_df, sat_annot[['chr', 'pos', 'alt'] + tp53_assays], on=['chr', 'pos', 'alt'], how='left')

            if gene == 'PTEN':
                cv_annot = pd.merge(cv_df, sat_annot[['chr', 'pos', 'alt'] + pten_assays], on=['chr', 'pos', 'alt'], how='left')
            
            if gene == 'KRAS':
                cv_annot = pd.merge(cv_df, sat_annot[['chr', 'pos', 'alt'] + ras_assays], on=['chr', 'pos', 'alt'], how='left')
            
            if gene in ['HRAS', 'NRAS']:
                cv_annot = pd.merge(cv_df, sat_annot[['chr', 'pos', 'alt'] + ['Ras_Bandaru']], on=['chr', 'pos', 'alt'], how='left')
            
            if gene == 'DNMT3A':
                cv_annot = pd.merge(cv_df, sat_annot[['chr', 'pos', 'alt'] + dnmt3a_assays], on=['chr', 'pos', 'alt'], how='left')
        
        else:

            cv_annot = cv_df.copy()

        # annotate dbNSFP data
            
        sat_annot_fn = os.path.join(os.environ['OUTPUT'], 'benchmarks', 'saturation_dbNSFP', f'{gene}.{ttype}.saturation.dbNSFP.tsv')
        sat_annot = pd.read_csv(sat_annot_fn, sep='\t', low_memory=False)
        sat_annot['chr'] = sat_annot['chr'].astype(str)

        cv_annot = pd.merge(cv_annot, sat_annot[['chr', 'pos', 'alt'] + all_dbnsfp_scores], on=['chr', 'pos', 'alt'], how='left')
        cv_annot.to_csv(f'{gene}.{ttype}.50.iter.annotated.tsv', sep='\t', index=False)


if __name__ == '__main__':

    cli()
