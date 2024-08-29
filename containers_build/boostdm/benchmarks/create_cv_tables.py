import click
import glob
import os
import tqdm
import gzip
import pickle

import pandas as pd

from boostdm import BoostDMError

# paths

basepath = os.environ['OUTPUT']
cvsplit_folder = os.path.join(basepath, 'splitcv_meta')  # cross-validation split folders
models_folder = os.path.join(basepath, 'training_meta')  # models folder
modeleval_folder = os.path.join(basepath, 'evaluation')  # cross-validation scores after training


def parse_cv(gene, ttype):
    
    fn = os.path.join(cvsplit_folder, ttype, f'{gene}.cvdata.pickle.gz')
    with gzip.open(fn, 'rb') as f:
        l = pickle.load(f)
    return l


def load_models(gene, ttype):
    
    path = os.path.join(models_folder, f'{ttype}/{gene}.models.pickle.gz')
    with gzip.open(path, 'rb') as g:
        models = pickle.load(g)
    return models


def get_cv_prediction(gene, ttype, exclude_train=False):
        
    models = load_models(gene, ttype)
    splitcv = parse_cv(gene, ttype)
    
    res = []

    for i, split in enumerate(splitcv):
        
        train_features = split[0]
        test_features = split[1]
        
        if exclude_train:
            # mutations in the test set, excluding those that are already seen in the training set
            muts = test_features[
                ~test_features.apply(
                    lambda x: (x['pos'], x['ref'], x['alt']), 
                    axis=1
                ).isin(train_features.apply(lambda x: (x['pos'], x['ref'], x['alt']), axis=1))
            ]
        else:
            muts = test_features.copy()
                    
        features = muts[[c for c in muts.columns if c not in ['chr', 'pos', 'ref', 'alt']]]
        coords = muts[['chr', 'pos', 'ref', 'alt']]
        
        y = split[3].loc[muts.index]  # true driver labels
        
        mod = models['models'][i]

        # prediction of i-th base model against muts
        yhat = mod.predict_proba(features)[:, 1]

        # put together features, genomic coordinates and true driver labels
        df = pd.concat([features, coords, y], axis=1)
                
        # add prediction with i-th base model
        df['boostDM_score'] = yhat
        df.rename(columns={'label': 'driver'}, inplace=True)
        
        res.append(df)

    return res


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
    
    visited = []

    for gene, ttype_model, _ in tqdm.tqdm(pool):

        if (gene, ttype_model) not in visited:

            res = get_cv_prediction(gene, ttype_model)
            
            # create table with the first iteration
            full_table = res[0]
            full_table = full_table.assign(iteration=0)

            for i in range(1, 50):
                iteration = res[i].assign(iteration=i)
                full_table = pd.concat([full_table, iteration])
            
            full_table = full_table.reset_index(drop=True)
            
            # save
            fn =  f'{gene}.{ttype_model}.50.iter.tsv'
            full_table.to_csv(fn, sep='\t', index=False)

            visited.append((gene, ttype_model))


if __name__ == '__main__':

    cli()
