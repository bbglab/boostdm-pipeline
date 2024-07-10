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


saturation_folder = os.path.join(os.environ['OUTPUT'], 'saturation', 'prediction')

cohorts = pd.read_csv(os.path.join(os.environ['INTOGEN_DATASETS'], 'cohorts.tsv'), sep='\t')
cohort_ttypes = cohorts['CANCER_TYPE'].unique()


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


# clustering options


def config_params(font_size=10):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


def round_low(x):
    
    y = abs(x - np.round(x, 3))
    return x - y


def plot_cluster_domains_kde(df, df_pfam, df_names, gene, output='./', dpi=150, invisible_heatmap=False, plot=True):
    
    X = df.values
    Y = pdist(X, metric=mcc_dist)

    if len(X) > 1:
        linkage = hierarchy.linkage(Y, method='ward')
    
    
    if gene == 'TP53':
        figsize = (13, 8)
    else:
        figsize = (13, 6)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,1], height_ratios=[1,16])
    gs.update(hspace=0.05, wspace=0.00)
    ax0 = plt.subplot(gs[0]) # counts_muts
    ax1 = plt.subplot(gs[1]) # null
    ax2 = plt.subplot(gs[2]) # heatmap
    ax3 = plt.subplot(gs[3]) # right dendogram
        
    # plot dendrogram and display cophenetic distances in ax3
    
    ax3.axis('off')
    
    if len(X) > 1:
        
        ddgram = dendrogram(linkage, truncate_mode=None,
                            labels=df.index,
                            color_threshold=0,
                            above_threshold_color='black',
                            orientation="right",
                            get_leaves=True,
                            no_plot=False, ax=ax3)
        
        # capture spread of y-values of dendrogram to adjust later
        
        pool = []
        for item in ddgram['icoord']:
            pool += item
        min_ = min(pool)
        max_ = max(pool)
        
        coph_diam = round_low(ddgram['dcoord'][-1][1])
    
    # Draw domains

    ax1.axis('off')
    ax0.axis('off')
    
    ax0.axhline(y=0.0, xmin=0, xmax=df.shape[1], ls="-", lw=2,color="black", alpha=0.5, zorder=1)
    
    fontsize = 7

    if df_pfam[df_pfam['SYMBOL'] == gene].shape[0] > 0:
        
        dg_pfam = get_PFAMs_per_transcript(gene, df_pfam, df_names)

        for i, r in dg_pfam.sort_values(by='START').iterrows():
            
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

            ax0.annotate(r["DOMAIN_NAME"], xy=(start_base + 2, 1), fontsize=fontsize, zorder=10)

            ax0.set_xlim(0, df.shape[1] + 50)
            ax0.set_ylim(-10, 10)
    
    # Heatmap
    
    y_values, x_values = [], []
    
    if len(X) > 1:
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
    
    eps = 0.05
    ax2.scatter(x_values, [y + np.random.uniform(-eps, eps) for y in y_values], color="#cc0000", s=1., alpha=0.5)

    # separating horizontal lines
    
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xlabel("Amino acid position",fontsize=12)
    ax2.set_ylabel("Tumor type",fontsize=12, rotation=90)
    
    if len(X) > 1:
        # labels = [str(ttype) for ttype in ddgram["ivl"][::-1]]
        labels = [str(ttype) for ttype in ddgram["ivl"]]

    else:
        labels = [str(ttype) for ttype in df.index]
    
    # len_cds = (df.columns.values[-1] + 1)
    len_cds = df.shape[1]
    len_aa =  len_cds

    ax2.set_xticks([x for x in np.linspace(0, len_cds, num=10, endpoint=True)])
    ax2.set_xticklabels([str(int(x)) for x in np.linspace(0, len_aa, num=10, endpoint=True)])
    
    if len(X) > 1:
        _ = ax2.set_yticks(np.linspace(min_, stop=max_, num=len(labels)))
        _ = ax2.set_yticklabels(labels, rotation=0, fontsize=10)
        h = ax2.hlines(np.linspace(min_, stop=max_, num=len(labels)), xmin=0, xmax=len_cds, alpha=0.3)
    else:
        _ = ax2.set_yticks(np.linspace(0, stop=len(X)+1, num=len(labels)))
        _ = ax2.set_yticklabels(labels, rotation=0, fontsize=10)
        h = ax2.hlines(np.linspace(0, stop=len(X)+1, num=len(labels)), xmin=0, xmax=len_cds, alpha=0.3)
    h.set_linewidth(0.5)
    
    ax2.set_xlim(0, len_cds)

    if gene == "TP53":
        ax2.tick_params(axis='y', labelsize=8.5, pad=0.25, width=0.5, length=1.5)
    else:
        ax2.tick_params(axis='both', labelsize=10, pad=0.25, width=0.5, length=1.5)
    
    if invisible_heatmap:
        ax2.set_visible(False)
     
    title = f'{gene}'
    ax0.set_title(title, fontsize=14)
    
    if len(X) > 1:
        ax3.set_ylim(ax2.get_ylim())
    else:
        ax2.set_ylim(-0.5, 0.5)
    
    plt.savefig(f'{gene}.clustered_blueprint.png', dpi=dpi, bbox_inches='tight')
    plt.savefig(f'{gene}.clustered_blueprint.svg', dpi=dpi, bbox_inches='tight')

    if not plot:
        plt.close(fig)
    

def create_saturation_table_reformat(pair_vectors):
    
    """
    pair_vectors: type (ttype, gene): array-like
    all genes are supposed to be the same
    """
    
    l_data = []
    l_ttype = []
    l_model = []
    
    for (ttype, gene), v in pair_vectors.items():
        
        l_data.append(list(v))
        l_ttype.append(ttype)

    df = pd.DataFrame(l_data)
    df.fillna(0.0, inplace=True)
    df.index = l_ttype
    return df


def plot_clustering_reformat(gene, df_pfam, df_names, output, res, plot=False, dpi=150, invisible_heatmap=False):
    
    """res: type (ttype, gene): (array-like, array-like)"""
    
    pair_vectors = {}
    
    for tt, g in res:
        if (g == gene) and (tt in cohort_ttypes):
            pair_vectors[(tt, gene)] = res[(tt, gene)]
    
    if len(pair_vectors) == 0:
        raise ValueError(f'{gene} does not have any blueprints')
    
    if len(pair_vectors) > 0:
        df = create_saturation_table_reformat(pair_vectors)
        x = plot_cluster_domains_kde(df, df_pfam, df_names, gene, plot=plot, output=output, dpi=dpi, invisible_heatmap=invisible_heatmap)



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
            plot_clustering_reformat(gene, df_pfam, df_names, './', saturation_vectors, dpi=300, plot=False)
        except ValueError as e:
            print(e)
    
    return


if __name__ == '__main__':
    
    cli()
