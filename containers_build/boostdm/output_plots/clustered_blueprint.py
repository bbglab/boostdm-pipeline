import sys
import os
import glob
import gzip
import tqdm
import pickle
import json
import click

import numpy as np
import pandas as pd
import scipy
import re

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from matplotlib import gridspec

from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import fcluster, cophenet
import scipy.cluster.hierarchy as hierarchy
from scipy.cluster.hierarchy import dendrogram

from boostdm import BoostDMError


saturation_folder = os.path.join(os.environ['OUTPUT'], 'saturation', 'prediction')


def get_aa_position(aachange):
    
    digit = re.findall(r'\d+', aachange)[0]
    return int(digit)


def create_saturation_vectors():

    res = {}

    for fn in tqdm.tqdm(glob.glob(os.path.join(saturation_folder, '*.prediction.tsv.gz'))):
        
        gene = os.path.basename(fn).split('.')[0]
        tt_model = os.path.basename(fn).split('.')[2]
        tt_features = os.path.basename(fn).split('.')[4]

        if tt_model != tt_features:
            continue

        df = pd.read_csv(fn, sep='\t')
        df = df[df['aachange'].notna()]
        df['aa_position'] = df['aachange'].apply(get_aa_position)
        df.sort_values(by=['aa_position'], inplace=True)
        dg = df.groupby('aa_position').agg({'boostDM_class': 'max'}).reset_index()

        saturation_vector = np.array(list(map(int, dg['boostDM_class'].values)))
        res[(tt_model, gene)] = saturation_vector
    
    return res


def get_PFAMs_per_transcript(gene, df_pfam, df_names):

    df_pfam_gene = df_pfam[(df_pfam['SYMBOL'] == gene)]
    df_pfam_gene = df_pfam_gene[['SYMBOL', 'START', 'END', 'DOMAIN']].drop_duplicates()
    df_pfam_gene = pd.merge(df_pfam_gene, df_names[['DOMAIN', 'DOMAIN_NAME']].drop_duplicates(), how="left")
    df_pfam_gene['POS'] = df_pfam_gene.apply(lambda row: row['START'] + (row['END'] - row['START']) // 2, axis=1)
    df_pfam_gene['SIZE'] = df_pfam_gene.apply(lambda row: row['END'] - row['START'] + 1, axis=1)
    df_pfam_gene['Color'] = "#998ec3"
    
    return df_pfam_gene


# distance used for clustering

def mcc_score(x, y):
    
    """Generalization of MCC score to [0, 1] continuous values"""
    
    x = np.array(x)
    y = np.array(y)
    
    tp = np.dot(x, y)
    tn = np.dot(1 - x, 1 - y)
    fp = np.dot(1 - x, y)
    fn = np.dot(x, 1 - y)

    # MCC = TP * TN – FP * FN / √ (TP +FP) * (TP + FN) * (TN + FP) * (TN + FN)
    
    num = (tp * tn) - (fp * fn) 
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    
    return num / den


def mcc_dist(x, y):
    
    return max(1 - mcc_score(x, y), 0)


def clustered_blueprints(gene, saturation_vectors, df_pfam, df_names):

    # create dataframe object with rows binary vector of aa positions indexed by tumor type

    subdict = {}
    for tt, g in saturation_vectors:
        if g == gene:
            subdict[(tt, gene)] = saturation_vectors[(tt, gene)]
    
    if len(subdict) == 0:
        raise BoostDMError(f'No blueprint for gene={gene}')

    l_data, l_ttype = [], []
    for (ttype, gene), v in subdict.items():        
        l_data.append(list(v))
        l_ttype.append(ttype)

    df = pd.DataFrame(l_data)
    df.fillna(0.0, inplace=True)
    df.index = l_ttype

    # define the main data: X driver mutations, Y distance matrix

    X = df.values
    Y = pdist(X, metric=mcc_dist)

    # define the plot environment with 4 axes
    
    if X.shape[0] > 20:

        fig_height = 6 + int(0.1 * X.shape[0])    
        fig = plt.figure(figsize=(13, fig_height))
        gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,1], height_ratios=[1,20])
        gs.update(hspace=0.05, wspace=0.00)

    elif X.shape[0] > 4:
        
        fig_height = 6 + int(0.1 * X.shape[0])    
        fig = plt.figure(figsize=(13, fig_height))
        gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,1], height_ratios=[1,13])
        gs.update(hspace=0.05, wspace=0.00)
    
    elif X.shape[0] > 1:

        fig_height = 3 + int(0.1 * X.shape[0])    
        fig = plt.figure(figsize=(13, fig_height))
        gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,1], height_ratios=[1,8])
        gs.update(hspace=0.05, wspace=0.00)

    else:

        fig_height = 2 + int(0.1 * X.shape[0])    
        fig = plt.figure(figsize=(13, fig_height))
        gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,1], height_ratios=[1,4])
        gs.update(hspace=0.05, wspace=0.00)

    
    ax0 = fig.add_subplot(gs[0])  # top-left: axis for domains
    ax1 = fig.add_subplot(gs[1])  # top-right: not used
    ax2 = fig.add_subplot(gs[2])  # bottom-left: blueprint tracks
    ax3 = fig.add_subplot(gs[3])  # bottom-right: dendrogram

    # perform hierarchical clustering with binary dataframe using MCC as distance

    X = df.values
    Y = pdist(X, metric=mcc_dist)
    
    if X.shape[0] > 1:

        linkage = hierarchy.linkage(Y, method='ward')
        ddgram = dendrogram(linkage, truncate_mode=None,
                            labels=df.index,
                            color_threshold=0,
                            above_threshold_color='black',
                            orientation="right",
                            get_leaves=True,
                            no_plot=False, ax=ax3)
                        
    ax3.axis('off')  # removes the names of the leaves in the dendrogram

    # draw Pfam domains

    ax1.axis('off')
    ax0.axis('off')
    
    ax0.axhline(y=0.0, xmin=0, xmax=df.shape[1], ls="-", lw=2,color="black", alpha=0.5, zorder=1)
    
    if df_pfam[df_pfam['SYMBOL'] == gene].shape[0] > 0:
        
        dg_pfam = get_PFAMs_per_transcript(gene, df_pfam, df_names)

        for i, r in dg_pfam.sort_values(by='START').reset_index().iterrows():
            
            start_base = r['START']
            size_base = r['SIZE']

            rect1 = patches.Rectangle(xy=(start_base, -1), width=size_base, 
                                        height=10, color="#90EE90", alpha=1, 
                                        clip_on=True, zorder=10)
                                        
            rect2 = patches.Rectangle(xy=(start_base, -1), width=size_base, 
                                        height=10, color="black", fill=None, alpha=1, 
                                        clip_on=True, zorder=10, linewidth=0.5)
            ax0.add_patch(rect1)
            ax0.add_patch(rect2)

            # write Pfam domain labels

            if gene in ['ARID2', 'ARHGAP35', 'AXIN1', 'BCL9L', 'BRCA1', 'CDH1', 'CTCF', \
                        'EZH2', 'FUBP1', 'KDM6A', 'PBRM1', 'RBM10', 'RNF43', 'SMARCA4', \
                        'WT1',  'XPO1']:
                y_min, y_max = -30, 70
                gap = y_max - y_min
                levels = [y_min, 
                          y_min + 0.4 * gap, 
                          y_min + 0.7 * gap, 
                          y_min + 1 * gap]
                ax0.annotate(r["DOMAIN_NAME"], xy=(start_base + 2, levels[i % 4]), fontsize=7, zorder=10, rotation=0)
                ax0.set_ylim(y_min, y_max)
            
            elif X.shape[1] < 2000:
                ax0.annotate(r["DOMAIN_NAME"], xy=(start_base + 2, 1), fontsize=7, zorder=10)
                ax0.set_ylim(-10, 10)
            else:
                y_min, y_max = -30, 70
                gap = y_max - y_min
                levels = [y_min, 
                          y_min + 0.4 * gap, 
                          y_min + 0.7 * gap, 
                          y_min + 1 * gap]
                ax0.annotate(r["DOMAIN_NAME"], xy=(start_base + 2, levels[i % 4]), fontsize=7, zorder=10, rotation=0)
                ax0.set_ylim(y_min, y_max)

            ax0.set_xlim(0, df.shape[1] + 50)

    # draw blueprints

    # gets the spread of the dendrogram

    if X.shape[0] > 1:
        pool = []
        for item in ddgram['icoord']:
            pool += item
        min_, max_ = min(pool), max(pool)
    else:
        min_, max_ = None, None
    
    y_values, x_values = [], []
    
    if X.shape[0] > 1:
        values = df.loc[ddgram["ivl"][::-1]].values
        scaling_factor = (max_ - min_) / (len(X) - 1)
    else:
        values = [df.iloc[0].values]
        scaling_factor = 1
        max_ = 0

    j = 0
    for array in values:
        for v in range(0, len(array)):
            if array[v] > 0:
                y_values.append(max_ + j * scaling_factor)
                x_values.append(v)
        j-=1

    # add jitter to y_values

    if X.shape[0] > 6:
        sigma_jitter = 0.75
    elif X.shape[0] > 1:
        sigma_jitter = 0.15
    else:
        sigma_jitter = 0.025

    y_values = np.array(y_values) + np.random.normal(0, sigma_jitter, size=len(y_values))

    ax2.scatter(x_values, y_values, color="#cc0000", s=1., alpha=0.5)
    
    # aminoacid positions

    len_aa = df.shape[1]
    ax2.set_xticks([x for x in np.linspace(0, len_aa, num=10, endpoint=True)])
    ax2.set_xticklabels([str(int(x)) for x in np.linspace(0, len_aa, num=10, endpoint=True)])    
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xlabel("Amino acid position",fontsize=12)
    ax2.set_ylabel("Tumor type",fontsize=12, rotation=90)

    # tumor type labels

    if X.shape[0] > 1:
        labels = [str(ttype) for ttype in ddgram["ivl"]]
    else:
        labels = [str(ttype) for ttype in df.index]


    if X.shape[0] > 1:
        _ = ax2.set_yticks(np.linspace(min_, stop=max_, num=len(labels)))
        _ = ax2.set_yticklabels(labels, rotation=0, fontsize=7)
        h = ax2.hlines(np.linspace(min_, stop=max_, num=len(labels)), xmin=0, xmax=len_aa, alpha=0.3)
    else:
        _ = ax2.set_yticks(np.linspace(0, stop=len(X)+1, num=len(labels)))
        _ = ax2.set_yticklabels(labels, rotation=0, fontsize=7)
        h = ax2.hlines(np.linspace(0, stop=len(X)+1, num=len(labels)), xmin=0, xmax=len_aa, alpha=0.3)
    
    h.set_linewidth(0.5)

    # embellishments

    ax2.set_xlim(0, len_aa)
    ax2.tick_params(axis='both', labelsize=10, pad=0.25, width=0.5, length=1.5)

    if X.shape[0] > 1:
        ax3.set_ylim(ax2.get_ylim())
    else:
        ax2.set_ylim(-0.5, 0.5)

    fig.suptitle(gene, fontsize=14)

    plt.savefig(f'{gene}.clustered_blueprint.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{gene}.clustered_blueprint.svg', dpi=300, bbox_inches='tight')
    plt.close()


@click.command()
def cli():

    saturation_vectors = create_saturation_vectors()

    PFAM_file = os.path.join(os.environ['BOOSTDM_DATASETS'], 'pfam_biomart.tsv.gz')
    PFAM_info = os.path.join(os.environ['BOOSTDM_DATASETS'], 'pfam_info.name.tsv')
    regions_fn = os.path.join(os.environ['BOOSTDM_DATASETS'], 'saturation' , 'cds-5spli.regions.gz')

    regions = pd.read_csv(regions_fn, sep='\t', usecols=[5, 6], low_memory=False).drop_duplicates()
    df_pfam = pd.read_csv(PFAM_file, sep="\t", names=["ENSEMBL_GENE", "TRANSCRIPT_ID", "START", "END", "DOMAIN"])
    df_names = pd.read_csv(PFAM_info, sep="\t", names=["DOMAIN", "CLAN", "CLAN_NAME", "DOMAIN_NAME", "Long Name"])

    df_pfam = df_pfam.merge(regions)
    df_pfam.head()  

    drivers = pd.read_csv(os.path.join(os.environ['INTOGEN_DATASETS'], 'drivers.tsv'), sep='\t')
    driver_genes = drivers['SYMBOL'].unique()

    for gene in tqdm.tqdm(driver_genes):
        try:
            clustered_blueprints(gene, saturation_vectors, df_pfam, df_names)
        except BoostDMError as e:
            print(e)
    
    return


if __name__ == '__main__':
    
    cli()
