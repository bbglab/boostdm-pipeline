import os
import glob
import tqdm
import pandas as pd

# global paths

boostdm_output = os.environ['OUTPUT']
saturation_mutagenesis = os.path.join(boostdm_output, 'saturation', 'annotation')


def vep_format(df):

    dg = pd.DataFrame()
    dg['chr'] = df['chr'].values
    dg['start'] = df['pos'].values
    dg['end'] = df['pos'].values
    dg['mut'] = df.apply(lambda r: f'{r["Reference"]}/{r["alt"]}', axis=1)
    dg = dg.sort_values(['chr', 'start'])
    
    return dg


if __name__ == '__main__':

    visited = []
    for fn in tqdm.tqdm(glob.glob(os.path.join(saturation_mutagenesis, '*.tsv.gz'))):
        gene = os.path.basename(fn).split('.')[0]
        if gene not in visited:
            df = pd.read_csv(fn, sep='\t', low_memory=False)
            dg = vep_format(df)
            dg.to_csv(f'{gene}.tsv', sep='\t', index=None, header=False)
            visited.append(gene)
