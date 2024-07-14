import os
import tqdm
import glob
import click

import pandas as pd
import numpy as np

from sklearn.metrics import precision_recall_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc

import matplotlib.pyplot as plt


# basepath

basepath = os.environ['OUTPUT']

# dbNSFP scores

all_dbnsfp_scores = ['AlphaMissense_score',
 'CADD_raw',
 'ESM1b_score',
 'EVE_score',
 'FATHMM_score',
 'MetaLR_score',
 'MetaRNN_score',
 'MutationAssessor_score',
 'PROVEAN_score',
 'Polyphen2_HDIV_score',
 'Polyphen2_HVAR_score',
 'REVEL_score',
 'SIFT4G_score',
 'SIFT_score',
 'VEST4_score',
 'phyloP100way_vertebrate',
 'phyloP17way_primate',
 'phyloP470way_mammalian']

# MAVE scores

tp53_assays = ['TP53_Boettcher', 'TP53_Giacomelli', 'TP53_Kato', 'TP53_Kotler', 'TP53_Ursu']
dnmt3a_assays = ['DNMT3A_Lue']
ras_assays = ['Ras_Bandaru', 'KRAS_Ursu']
pten_assays = ['PTEN_Mighell', 'PTEN_Matreyek']

# plotting palette

palette_boostdm = {'boostDM_score': '#ac0f0f'}

palette_dnNSFP = {'AlphaMissense_score': '#FFAD84', 'EVE_score': '#FFC47E', 'ESM1b_score': '#FFE382', 
                  'CADD_raw': '#E26EE5', 'FATHMM_score': '#6C22A6', 'MetaLR_score': '#6962AD', 'MetaRNN_score': '#A084E8',
                  'MutationAssessor_score': '#EA8FEA', 'PROVEAN_score': '#FFAACF', 
                  'Polyphen2_HDIV_score': '#C21292', 'Polyphen2_HVAR_score': '#711DB0', 
                  'REVEL_score': '#5D3587', 'SIFT4G_score': '#FED9ED', 'SIFT_score': '#E7BCDE', 'VEST4_score': '#BB9CC0',
                  'phyloP100way_vertebrate': '#944E63', 'phyloP17way_primate': '#B47B84', 'phyloP470way_mammalian': '#CAA6A6'
                  }

palette_mave = {'TP53_Kato': '#294B29', 
                'TP53_Giacomelli': '#50623A', 
                'TP53_Kotler': '#789461', 
                'TP53_Ursu': '#DBE7C9',
                'TP53_Boettcher': '#E1F0DA',
                
                'PTEN_Mighell': '#83C0C1',
                'PTEN_Matreyek': '#96E9C6',

                'Ras_Bandaru': '#43766C',
                'KRAS_Ursu': '#607274',
                
                'DNMT3A_Lue': '#113946'}

palette = {**palette_boostdm, **palette_dnNSFP, **palette_mave}

phylop_scores = ['phyloP100way_vertebrate', 'phyloP17way_primate', 'phyloP470way_mammalian']
generative_scores = ['AlphaMissense_score', 'EVE_score', 'ESM1b_score']
classical_scores = ['CADD_raw', 'FATHMM_score', 'MetaLR_score', 'MetaRNN_score',
                  'MutationAssessor_score', 'PROVEAN_score', 
                  'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 
                  'REVEL_score', 'SIFT4G_score', 'SIFT_score', 'VEST4_score']

score_dict = {
    'phylop': phylop_scores,
    'generative': generative_scores,
    'classical': classical_scores,
    'tp53': tp53_assays,
    'dnmt3a': dnmt3a_assays,
    'pten': pten_assays,
    'ras': ras_assays
}


def plot_prc(gene, ttype, score, ax=None, plot=True, **kwargs):

    cv_table = os.path.join(basepath, 'benchmarks', 'cv_tables_annotated', f'{gene}.{ttype}.50.iter.annotated.tsv')

    try:
        df = pd.read_csv(cv_table, sep='\t')
    except Exception as e:
        print(gene, ttype, e)


    if not score in df.columns:
        return
    
    df = df[['driver', score]].dropna()

    # make the mutation set balanced

    prop_drivers = df['driver'].mean()

    if prop_drivers >= 0.5:
        num_to_drop = df[df['driver'] == 1].shape[0] - df[df['driver'] == 0].shape[0]
        if num_to_drop > 0:
            drop_indices = np.random.choice(df[df['driver'] == 1].index, num_to_drop, replace=False)
            df = df.drop(drop_indices)
    else:
        num_to_drop = df[df['driver'] == 0].shape[0] - df[df['driver'] == 1].shape[0]
        if num_to_drop > 0:
            drop_indices = np.random.choice(df[df['driver'] == 0].index, num_to_drop, replace=False)
            df = df.drop(drop_indices)
    
    if df.shape[0] == 0:
        return

    # fit simple logistic model

    model = LogisticRegression(solver='lbfgs')
    X = df[score].values.reshape(-1, 1)
    y = df['driver'].values
    
    # number of mutations
    
    n = len(y)
    positive = sum(y)
    negative = n - positive
    
    model.fit(X, y)
    yhat = model.predict_proba(X)
    probs = yhat[:, 1]
    assert(X.shape[0] == len(y))

    # precision-recall auc
    
    precision, recall, _ = precision_recall_curve(y, probs)
    auc_score = auc(recall, precision)

    # compute and return precision-recall curves
    # only if there is a minimum number of positive instances

    if positive > 20:
        if plot:
            plot_pr_curve(y, probs, ax, label=f'{score}: auPRC={auc_score:.2}', **kwargs)
        return auc_score, precision, recall, positive, negative
    
    return


def plot_pr_curve(testy, model_probs, ax, **kwargs):
    """plot model precision-recall curve"""

    precision, recall, _ = precision_recall_curve(testy, model_probs)
    ax.plot(recall, precision, **kwargs)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_ylim(0.3, 1.01)
    ax.set_xlim(0.01, 1.01)


def plot_several(gene, ttype, scores, save_path=None):
    
    fig, ax = plt.subplots(figsize=(4,4))
    
    for s in scores:
        plot_prc(gene, ttype, s, ax=ax, color=palette[s], lw=5, alpha=0.5)
    plot_prc(gene, ttype, 'boostDM_score', ax=ax, color=palette['boostDM_score'], lw=5, alpha=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(f'{gene} ({ttype})')
    ax.legend(loc=(0.5, 0.01))
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.show()

def plot_score_family(gene, ttype, score_family):
    
    fig, ax = plt.subplots(figsize=(4,4))

    all_none = True
    for s in score_dict[score_family]:
        arg = plot_prc(gene, ttype, s, ax=ax, color=palette[s], lw=5, alpha=0.5)
        if arg is not None:
            all_none = False
    
    if all_none:
        return

    plot_prc(gene, ttype, 'boostDM_score', ax=ax, color=palette['boostDM_score'], lw=5, alpha=1)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(f'{gene} ({ttype})')
    ax.legend(loc=(0.5, 0.01))
    
    fig.savefig(f'{gene}.{ttype}.{score_family}.prauc.svg', bbox_inches='tight')
    fig.savefig(f'{gene}.{ttype}.{score_family}.prauc.png', dpi=150, bbox_inches='tight')
    
    plt.close()


@click.command()
@click.option("--plots", is_flag=True, show_default=True, default=False, help="Dump all PR plots")
def cli(plots):

    pool = []
    for fn in glob.glob(os.path.join(basepath, 'benchmarks', 'cv_tables_annotated', f'*.*.50.iter.annotated.tsv')):
        gene = os.path.basename(fn).split('.')[0]
        ttype = os.path.basename(fn).split('.')[1]
        pool.append((gene, ttype))
    
    if plots:

        # SCORE FAMILIES
        # dump all PR plots

        for gene, ttype in tqdm.tqdm(pool):
            for score_family in score_dict:
                plot_score_family(gene, ttype, score_family)
    
    else:

        # ALL SCORES
        # dump all PR-AUC in a table


        scores = all_dbnsfp_scores + tp53_assays + dnmt3a_assays + ras_assays + pten_assays + ['boostDM_score']
        cols = ['gene', 'ttype'] + scores
        res_dict = {}

        for gene, ttype_model in tqdm.tqdm(pool):
            
            res_dict['gene'] = res_dict.get('gene', []) + [gene]
            res_dict['ttype'] = res_dict.get('ttype', []) + [ttype_model]
            
            for s in scores:

                out = plot_prc(gene, ttype_model, s, ax=None, plot=False)
                if out is None:
                    res_dict[s] = res_dict.get(s, []) + [None]
                else:
                    auc_score, *_ = out
                    res_dict[s] = res_dict.get(s, []) + [auc_score]
        
        res_df = pd.DataFrame(res_dict)
        res_df.to_csv('benchmark.tsv', sep='\t', index=False)


if __name__ == '__main__':

    cli()

    """

    # MAVE: TP53

    for fn in tqdm.tqdm(glob.glob(os.path.join(basepath, 'benchmarks', 'cv_tables_annotated', '*.*.50.iter.annotated.tsv'))):
        gene, ttype = os.path.basename(fn).split('.')[:2]
        if gene == 'TP53':
            plot_several(gene, ttype, ['TP53_Kato', 'TP53_Giacomelli', 'TP53_Kotler', 'TP53_Ursu'])

    # MAVE: PTEN
    
    for fn in tqdm.tqdm(glob.glob(os.path.join(basepath, 'benchmarks', 'cv_tables_annotated', '*.*.50.iter.annotated.tsv'))):
        gene, ttype = os.path.basename(fn).split('.')[:2]
        if gene == 'PTEN':
            plot_several(gene, ttype, ['PTEN_Matreyek', 'PTEN_Mighell'])


    # MAVE: [KRAS, HRAS, NRAS]
    
    for fn in tqdm.tqdm(glob.glob(os.path.join(basepath, 'benchmarks', 'cv_tables_annotated', '*.*.50.iter.annotated.tsv'))):
        gene, ttype = os.path.basename(fn).split('.')[:2]
        if gene in ['HRAS', 'NRAS']:
            plot_several(gene, ttype, ['RAS_Bandaru'])
        if gene == 'KRAS':
            plot_several(gene, ttype, ['RAS_Bandaru', 'KRAS_Ursu'])

    # MAVE: DNMT3A
    
    for fn in tqdm.tqdm(glob.glob(os.path.join(basepath, 'benchmarks', 'cv_tables_annotated', '*.*.50.iter.annotated.tsv'))):
        gene, ttype = os.path.basename(fn).split('.')[:2]
        if gene == 'DNMT3A':
            plot_several(gene, ttype, ['DNMT3A_Lue'])

    """

