import os
import click
import glob
import tqdm

import pandas as pd


# config to read through MAVE files

tp53_assays = ['TP53_Boettcher', 'TP53_Giacomelli', 'TP53_Kato', 'TP53_Kotler', 'TP53_Ursu']
dnmt3a_assays = ['DNMT3A_Lue']
ras_assays = ['Ras_Bandaru', 'KRAS_Ursu']
pten_assays = ['PTEN_Mighell', 'PTEN_Matreyek']

mave_code = {
    'TP53_Kato': {
        'scores': ['tp53_score'],
        'aachange': 'aachange'
        },
    'TP53_Boettcher': {
        'scores': ['FUNC_SCORE'],
        'aachange': 'aachange'
        },
    'TP53_Ursu': {
        'scores': ['HotellingT2'],
        'aachange': 'Variant'
    },
    'TP53_Kotler': {
        'scores': ['RFS_H1299'],
        'aachange': 'aachange'
    },
    'TP53_Giacomelli': {
        'scores': ['FUNC_SCORE'],
        'aachange': 'aachange'
    },
    'Ras_Bandaru': {
        'scores': ['FUNC_SCORE'],
        'aachange': 'aachange'
    },
    'DNMT3A_Lue': {
        'scores': ['sgRNA_score_d9_citrine_positive'],
        'aachange': 'aachange'
    },
    'PTEN_Matreyek': {
        'scores': ['score'],
        'aachange': 'aachange'
    },
    'PTEN_Mighell': {
        'scores': ['Cum_score'],
        'aachange': 'aachange'
    },
    'KRAS_Ursu': {
        'scores': ['HotellingT2'],
        'aachange': 'Variant'
    }
}


# output folder

basepath = os.environ['OUTPUT']
mave_data_folder = os.environ['MAVE_DATA']


# command line interactive

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
        sat_df.reset_index(inplace=True, drop=True)
        sat_df['chr'] = sat_df['chr'].astype(str)
        sat_df['pos'] = sat_df['pos'].astype(int)
        sat_annot_df = sat_df.copy()

        if gene == 'TP53':

            for assay in tp53_assays:

                fn = os.path.join(mave_data_folder, assay, 'scores.tsv')
                annot_df = pd.read_csv(fn, sep='\t')
                aachange = mave_code[assay]['aachange']
                score    = mave_code[assay]['scores'][0]
                sat_annot_df = pd.merge(sat_annot_df, annot_df[[aachange, score]], how='left', left_on=['aachange'], right_on=[aachange], 
                            suffixes=('', '_y'))
                sat_annot_df = sat_annot_df.rename(columns={score: assay})
                sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)

            sat_annot_df.to_csv(f'{gene}.{ttype_model}.saturation.mave.tsv', sep='\t', index=False)

        if gene == 'PTEN':
            
            for assay in pten_assays:

                fn = os.path.join(mave_data_folder, assay, 'scores.tsv')
                annot_df = pd.read_csv(fn, sep='\t')
                aachange = mave_code[assay]['aachange']
                score    = mave_code[assay]['scores'][0]
                sat_annot_df = pd.merge(sat_annot_df, annot_df[[aachange, score]], how='left', left_on=['aachange'], right_on=[aachange], 
                            suffixes=('', '_y'))
                sat_annot_df = sat_annot_df.rename(columns={score: assay})
                sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)

            sat_annot_df.to_csv(f'{gene}.{ttype_model}.saturation.mave.tsv', sep='\t', index=False)

        if gene == 'DNMT3A':
            
            for assay in dnmt3a_assays:

                fn = os.path.join(mave_data_folder, assay, 'scores.tsv')
                annot_df = pd.read_csv(fn, sep='\t')
                aachange = mave_code[assay]['aachange']
                score    = mave_code[assay]['scores'][0]
                sat_annot_df = pd.merge(sat_annot_df, annot_df[[aachange, score]], how='left', left_on=['aachange'], right_on=[aachange], 
                            suffixes=('', '_y'))
                sat_annot_df = sat_annot_df.rename(columns={score: assay})
                sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)
                
            sat_annot_df.to_csv(f'{gene}.{ttype_model}.saturation.mave.tsv', sep='\t', index=False)
        
        if gene == 'KRAS':
            
            # KRAS Ursu
            assay = 'KRAS_Ursu'
            fn = os.path.join(mave_data_folder, assay, 'scores.tsv')
            annot_df = pd.read_csv(fn, sep='\t')

            aachange = mave_code[assay]['aachange']
            score    = mave_code[assay]['scores'][0]
            
            sat_annot_df = pd.merge(sat_annot_df, annot_df[[aachange, score]], how='left', left_on=['aachange'], right_on=[aachange], 
                        suffixes=('', '_y'))
            sat_annot_df = sat_annot_df.rename(columns={score: assay})
            sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)

            # Ras Bandaru
            assay = 'Ras_Bandaru'
            fn = os.path.join(mave_data_folder, assay, 'scores.tsv')
            annot_df = pd.read_csv(fn, sep='\t')
            annot_df = annot_df[annot_df['SYMBOL'] == gene]
            
            aachange = mave_code[assay]['aachange']
            score    = mave_code[assay]['scores'][0]
            
            sat_annot_df = pd.merge(sat_annot_df, annot_df[[aachange, score]], how='left', left_on=['aachange'], right_on=[aachange], 
                        suffixes=('', '_y'))
            sat_annot_df = sat_annot_df.rename(columns={score: assay})
            sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)

            sat_annot_df.to_csv(f'{gene}.{ttype_model}.saturation.mave.tsv', sep='\t', index=False)
        
        if gene in ['HRAS', 'NRAS']:

            assay = 'Ras_Bandaru'
            fn = os.path.join(mave_data_folder, assay, 'scores.tsv')
            annot_df = pd.read_csv(fn, sep='\t')
            annot_df = annot_df[annot_df['SYMBOL'] == gene]

            aachange = mave_code[assay]['aachange']
            score    = mave_code[assay]['scores'][0]
            
            sat_annot_df = pd.merge(sat_annot_df, annot_df[[aachange, score]], how='left', left_on=['aachange'], right_on=[aachange], 
                        suffixes=('', '_y'))
            sat_annot_df = sat_annot_df.rename(columns={score: assay})
            sat_annot_df.drop(sat_annot_df.filter(regex='_y$').columns, axis=1, inplace=True)

            sat_annot_df.to_csv(f'{gene}.{ttype_model}.saturation.mave.tsv', sep='\t', index=False)

        else:

            print(f'No MAVE data for gene={gene}')

if __name__ == '__main__':
    
    cli()
