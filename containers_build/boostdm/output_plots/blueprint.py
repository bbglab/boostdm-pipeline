from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm
from matplotlib import collections as mc
import os
import gzip
import pickle
import click


cancer_predictions_path = os.path.join(os.environ['OUTPUT'], 'saturation', 'prediction')
obs_mut = pd.read_csv(os.path.join(os.environ['OUTPUT'], 'discovery', 'mutations.tsv.gz'), sep='\t')
oncotree_table = pd.read_csv(os.path.join(os.environ['BOOSTDM_DATASETS'], 'shared', 'tree.tsv'), sep='\t')
cohorts = pd.read_csv(os.path.join(os.environ['INTOGEN_DATASETS'], 'cohorts.tsv'), sep='\t')

cmap = cm.get_cmap('tab10')
colors = ['#3a5a40', cmap(2), cmap(1), cmap(6), cmap(4), cmap(5), '#0077b6', '#03045e', '#00b4d8', '#90e0ef']

alphas = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]

names = {'HotMaps': '3D cluster',
         'CLUSTL': 'Linear cluster',
         'smRegions': 'PFAM domain',
         'PhyloP': 'Conservation',
         'PTM': 'PTM',
         'nmd': 'NMD',
         'csqn_type_missense': 'Missense',
         'csqn_type_nonsense': 'Nonsense',
         'csqn_type_splicing': 'Splicing'}

def get_PFAMs_per_transcript(transcript):
    
    BOOSTDM_DATASETS = os.environ['BOOSTDM_DATASETS']
    PFAM_DOMAINS_FILE = os.path.join(BOOSTDM_DATASETS, 'pfam_biomart.tsv.gz')
    PFAM_DOMAINS_INFO = os.path.join(BOOSTDM_DATASETS, 'pfam_info.name.tsv')
    
    df_pfam = pd.read_csv(PFAM_DOMAINS_FILE, sep="\t", names=["ENSEMBL_GENE", "ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"])
    df_names = pd.read_csv(PFAM_DOMAINS_INFO, sep="\t", names=["DOMAIN", "CLAN", "CLAN_NAME", "DOMAIN_NAME", "Long Name"])

    # Get domains
    df_pfam_gene = df_pfam[(df_pfam["ENSEMBL_TRANSCRIPT"] == transcript)]
    df_pfam_gene = df_pfam_gene[["ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"]].drop_duplicates()
    df_pfam_gene = pd.merge(df_pfam_gene, df_names[["DOMAIN", "DOMAIN_NAME"]].drop_duplicates(), how="left")
    if df_pfam_gene.shape[0] > 0:
        df_pfam_gene["POS"] = df_pfam_gene.apply(lambda row: row["START"] + ((row["END"] - row["START"]) // 2), axis=1)
        df_pfam_gene["SIZE"] = df_pfam_gene.apply(lambda row: row["END"] - row["START"] + 1, axis=1)
        df_pfam_gene["Color"] = "#808080ff"  # color Pfam domain

    return df_pfam_gene


"""Fetch mutations"""


def get_position(row):
    try:
        v = int("".join(row["aachange"][1:-1]))
        return v
    except Exception:
        return -1
    
def get_all_children(df, parent_id, included=None):
    if included is None:
        included = []

    included.append(parent_id)
    children = df[df['PARENT'] == parent_id]['ID'].tolist()

    for child_id in children:
        get_all_children(df, child_id, included)

    return included

    
def load_saturation_cancer(path, gene, shap_corrected):

    df = path.copy()
    df.drop_duplicates(inplace=True)
    df["Protein_position"] = df.apply(lambda row: get_position(row), axis=1)

    # aggregate PTMs
    df["PTM"] = df.apply(lambda r: max(r['Phosphorylation'],
                                           r['Acetylation'],
                                           r['Methylation'],
                                           r['Ubiquitination'],
                                           r['Regulatory_Site']), axis=1)

    # summarize to codon level
    df = df[~df['aachange'].isnull()]

    if len(df) > 0:
        df["AA"] = df.apply(lambda row: row["aachange"][0], axis=1)
        df3 = df.groupby(["AA", "Protein_position", "ENSEMBL_TRANSCRIPT"],
                         as_index=False).agg(
                            {
                                "boostDM_score": list,
                                "boostDM_class": np.any,
                                "HotMaps": np.nanmax,
                                "smRegions": np.nanmax,
                                "CLUSTL": np.nanmax,
                                "csqn_type_missense": list,
                                "csqn_type_nonsense": list,
                                "csqn_type_splicing": list,
                                "PhyloP": np.nanmean,
                                "PTM": np.nanmax,
                                "nmd": np.nanmax
                            })
        df3["gene"] = gene                                        
        df3.sort_values(by='Protein_position', ascending=True, inplace=True)
        df3.reset_index(inplace=True)
        return df3
    print(f"file {path_file} not found...")
    return pd.DataFrame([])


def create_observed_dataset(prediction_path, gene, cohort, obs_mut):
    sat_pred = pd.read_csv(prediction_path, sep='\t')
    sat_pred['chr'] = sat_pred['chr'].astype(str)
    sat_pred['pos'] = sat_pred['pos'].astype(int)
    
    ttypes_included = get_all_children(oncotree_table, cohort)
    matching_cohorts = []
    for ttype in ttypes_included:
        matching_cohorts += cohorts[cohorts['CANCER_TYPE'] == ttype]['COHORT'].tolist()
    matching_cohorts = set(matching_cohorts)

    obs_mut = obs_mut[(obs_mut['gene'] == gene) & (obs_mut['COHORT'].isin(matching_cohorts))].reset_index(drop=True)

    df = obs_mut.merge(sat_pred, on=['gene', 'chr', 'pos', 'alt', 'aachange'], how='left')
    df.rename(columns={"COHORT": "cancer_type"}, inplace=True)
    df = df[~df['boostDM_class'].isnull()]
    return df


def get_plot_data_joanen(data):

    data["Protein_position"] = data.apply(lambda row: get_position(row), axis=1)
    data["AA"] = data.apply(lambda row: row["aachange"][0], axis=1)
    data['ID'] = data.apply(lambda x: '{}_{}'.format(x['pos'], x['alt']), axis=1)
    data = data.groupby(["ID", "pos", "AA", "Protein_position", "gene", "ENSEMBL_TRANSCRIPT", "boostDM_score", "boostDM_class"], as_index=False).agg({"sampleID": "count"})
    data.rename(columns={"sampleID": "number_observed_muts"}, inplace=True)

    return data

def plot_gene_full_nucleotide(data, transcript, sat_pred, ax0, all_possible=False):

    # remove those mutations not falling in CDS:
    df = data[data['AA'] != 'n']

    # Configure the axis
    ax0.set_title('Observed mutations (n='+str(data['number_observed_muts'].sum())+')',  fontsize=6)
    ax0.set_ylabel("mutation count",  fontsize=6,  labelpad=3)

    ax0.spines['bottom'].set_linewidth(1)
    ax0.spines['left'].set_linewidth(1)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.tick_params(axis='y', labelsize=6, pad=3.5, width=0.75, length=3.5)
    ax0.tick_params(axis='x', length=0)
    ax0.set_xticks([])

   
    # set equivalent coordinates for the three possible mutations
    prot_pos = list(df.Protein_position)
    ys = df["number_observed_muts"].values
    d = df["boostDM_score"].values

    coordinates_mutations = []

    passenger_x = []
    passenger_y = []
    passenger_color = []

    driver_x = []
    driver_y = []
    driver_color = []

    # for each of the positions
    for i, p in enumerate(prot_pos):
        if ys[i] > 0:

            coordinates_mutations.append([(p, 0), (p, ys[i] - 0.1)])

            if d[i] < 0.5:

                passenger_x.append(p)
                if all_possible:
                    passenger_y.append(d[i])
                else:
                    passenger_y.append(ys[i])
                passenger_color.append('#636363')

            else:
                driver_x.append(p)
                if all_possible:
                    driver_y.append(d[i])
                else:
                    driver_y.append(ys[i])
                driver_color.append('#ac0f0f')

    lc = mc.LineCollection(coordinates_mutations, colors='black', linewidths=1, alpha=0.3)
    ax0.add_collection(lc)

    size = 12
    ax0.scatter(passenger_x, passenger_y, s=size, c=passenger_color, alpha=0.7, label='non-driver')
    ax0.scatter(driver_x, driver_y, s=size, c=driver_color, alpha=0.7, label='driver')

    leg = ax0.legend(loc=(0, 1.15), prop=dict(size=6))
    leg.get_frame().set_linewidth(0.0)

    ax0.set_xlim(0, max(sat_pred['Protein_position']))


def plot_codon_bands(df_pfam_gene, df, ax_0, ax_2, ax_4):

    ax_0.set_ylabel("boostDM score", fontsize=6)
    ax_0.set_xticks(np.linspace(0, 1, 3))
    ax_0.spines['bottom'].set_visible(False)
    ax_0.spines['left'].set_linewidth(1)
    ax_0.spines['right'].set_visible(False)
    ax_0.spines['top'].set_visible(False)

    # set equivalent coordinates for the three possible mutations
    prot_pos = list(df.Protein_position)
    
    # retrieve values
    d = df["boostDM_score"].values
    csqn_type_miss  = df["csqn_type_missense"].values 
    csqn_type_non   = df["csqn_type_nonsense"].values

    passenger_x, passenger_y, passenger_color = [], [], []
    driver_x, driver_y, driver_color = [], [], []

    # for each of the positions
    for i, p in enumerate(prot_pos):
        for j, score in enumerate(d[i]):
            
            if score < 0.5:
                passenger_x.append(p)
                passenger_y.append(score)
                passenger_color.append('#636363')
            else:
                driver_x.append(p)
                driver_y.append(score)
                driver_color.append('#800000')

    size = 1
    ax_0.scatter(passenger_x, passenger_y, s=size+2, c=passenger_color, alpha=0.1)
    ax_0.scatter(driver_x, driver_y, s=size+2, c=driver_color, alpha=0.1)
    ax_0.set_xticks([])
    ax_0.set_xlim(0, len(prot_pos))
    
    driver_score = []
    for aa, bDMscore in zip(df['Protein_position'], df['boostDM_score']):
        driver_score = driver_score+len([1 for x in bDMscore if float(x)>0.5])*[aa]
    
    miss_high_score = []
    for aa, bDMscore, ttype in zip(df['Protein_position'], df['boostDM_score'], df['csqn_type_missense']):
        ttype_score = [float(a)*float(b) for a,b in zip(bDMscore,ttype)]
        miss_high_score = miss_high_score+len([1 for x in ttype_score if float(x)>=0.9])*[aa]
    
    miss_low_score = []
    for aa, bDMscore, ttype in zip(df['Protein_position'], df['boostDM_score'], df['csqn_type_missense']):
        ttype_score = [float(a)*float(b) for a,b in zip(bDMscore,ttype)]
        miss_low_score = miss_low_score+len([1 for x in ttype_score if float(x)>=0.5 and float(x)<0.9])*[aa]

    non_high_score = []
    for aa, bDMscore, ttype in zip(df['Protein_position'], df['boostDM_score'], df['csqn_type_nonsense']):
        ttype_score = [float(a)*float(b) for a,b in zip(bDMscore,ttype)]
        non_high_score = non_high_score+len([1 for x in ttype_score if float(x)>=0.9])*[aa]
        
    non_low_score = []
    for aa, bDMscore, ttype in zip(df['Protein_position'], df['boostDM_score'], df['csqn_type_nonsense']):
        ttype_score = [float(a)*float(b) for a,b in zip(bDMscore,ttype)]
        non_low_score = non_low_score+len([1 for x in ttype_score if float(x)>=0.5 and float(x)<0.9])*[aa]

    ax_2.plot(driver_score, np.full_like(driver_score,5), '|k', markeredgewidth=0.2, markersize=6, color='#800000')
    ax_2.plot(miss_high_score, np.full_like(miss_high_score,4), '|k', markeredgewidth=0.2, markersize=6, color='#C17E14')
    ax_2.plot(miss_low_score, np.full_like(miss_low_score,3), '|k', markeredgewidth=0.2, markersize=6, color='#C17E14')
    ax_2.plot(non_high_score, np.full_like(non_high_score,2), '|k', markeredgewidth=0.2,  markersize=6, color='#C17E14')
    ax_2.plot(non_low_score, np.full_like(non_low_score,1), '|k', markeredgewidth=0.2,  markersize=6, color='#C17E14')

    # First number is the number of tiers that you have
    ax_2.set_yticks([1, 2, 3, 4, 5])
    ax_2.set_yticklabels(['nonsense low', 'nonsense high', 'missense low', 'missense high', 'all drivers'], fontsize = 6)
    ax_2.set_ylim(.5, 5.5)
    ax_2.set_xlim(0, max(df['Protein_position']))
    ax_2.set_xticks([])
    ax_2.spines['bottom'].set_visible(False)
    ax_2.spines['left'].set_linewidth(1)
    ax_2.spines['right'].set_visible(False)
    ax_2.spines['top'].set_linewidth(1)

    ax_4.set_ylim(0, 1)
    d = df["boostDM_score"].values

    for i, r in df_pfam_gene.iterrows():
        start_base = r['START']
        size_base = r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax_4.annotate(r["DOMAIN_NAME"], xy=(start_base + size_base/3, 0.3), fontsize=5)
        ax_4.add_patch(rect)

    protein_ticks = [0] + list(map(int, np.linspace(0, max(prot_pos), num=10)))

    ax_4.set_xticks(protein_ticks)
    ax_4.set_xticklabels(protein_ticks, fontsize = 6)
    
    ax_4.set_xlim(0, max(prot_pos))
    ax_4.set_yticks([])
    ax_4.tick_params(axis='x', which='major', pad=3)

    ax_0.set_yticks([0, 0.5, 1])
    ax_0.set_yticklabels([0, 0.5, 1], size=6)


def tracked_blueprint_all(gene, ttype_model, ttype_features, df_codon, df, sat_pred, show=False):
    wanted_df = df_codon[(df_codon['gene'] == gene)]

    for transcript, gene in wanted_df[["ENSEMBL_TRANSCRIPT", "gene"]].drop_duplicates().values:

        # get PFAM domains and subset the mutation data
        subset_data_pfam = get_PFAMs_per_transcript(transcript)

        subset_data_muts = df_codon[
            (df_codon["ENSEMBL_TRANSCRIPT"] == transcript)].sort_values(by='Protein_position', ascending=True)

        # define figure layout
        fig = plt.figure(figsize=(8, 5), dpi=300)

        # grid layout
        gs = gridspec.GridSpec(32, 13, figure=fig)

        border = -(len(names) + 1)  # bottom limit for blueprint scatter

        # define subplots space
        ax0 = plt.subplot(gs[:border-14, :12])
        ax1 = plt.subplot(gs[border-14:border-8, :12])
        ax3 = plt.subplot(gs[border-8:border-2, :12])
        ax5 = plt.subplot(gs[border-2:border-1, :12])

        axes = []
        for i, track in enumerate(names):
            axes.append(plt.subplot(gs[border+i+1, :12], sharex=ax1))
        
        # plot scatterplot, driver tracks and protein body
        plot_codon_bands(subset_data_pfam, subset_data_muts, ax1, ax3, ax5)

        # plot needlplot
        data = get_plot_data_joanen(df) #data, count_driver, count_total
        plot_gene_full_nucleotide(data, transcript, sat_pred, ax0)
        
        # plot features tracks
        for i, track in enumerate(names):
            color = colors[i]
            alpha = alphas[i]
            if track not in ['csqn_type_missense', 'csqn_type_nonsense', 'csqn_type_splicing', 'csqn_type_synonymous']:
                axes[i].plot(subset_data_muts[track].values, color=color, alpha=alpha, lw=0.5)
            else:
                datavector = np.array(list(map(np.nanmean, subset_data_muts[track].values)))
                axes[i].plot(datavector, color=color, alpha=alpha, lw=0.5)

            axes[i].spines['bottom'].set_visible(False)
            axes[i].spines['left'].set_linewidth(1)
            axes[i].spines['right'].set_visible(False)
            axes[i].spines['top'].set_visible(False)
            axes[i].set_yticks([])
            axes[i].set_xticks([])

            axes[i].set_ylabel(names[track], rotation=0, labelpad=8, fontsize=6, color=color,
                               horizontalalignment='right', verticalalignment='center')

        fn = os.path.join(os.path.join(os.environ['OUTPUT'], 'evaluation', ttype_model, f'{gene}.eval.pickle.gz'))
        
        with gzip.open(fn, 'rb') as f:
            l = pickle.load(f)
        Fscore50 = round(np.nanmean(l['fscore50']), 2)

        df_DI = pd.read_csv(os.path.join(os.environ['OUTPUT'], 'discovery', 'discovery.tsv.gz'), sep='\t')
        DiscoveryI = round(df_DI[(df_DI['gene']==gene) & (df_DI['ttype']==ttype_model)]['discovery_index'].iloc[0], 2)

        ax0.set_title(
            f'{gene} ({ttype_model}) \n' \
            f'{ttype_model} model (F-score50={str(Fscore50)}, discovery={str(DiscoveryI)})\n' \
            f'{ttype_model} mutations and features (n={str(len(df))})\n', fontsize=10, y=1.25)

        fn_svg = os.path.join(f'{gene}.model.{ttype_model}.features.{ttype_features}.svg')
        fn_png = os.path.join(f'{gene}.model.{ttype_model}.features.{ttype_features}.png')
        plt.savefig(fn_svg, bbox_inches='tight')
        plt.savefig(fn_png, dpi=300, bbox_inches='tight')


@click.command()
@click.option('--gene', type=str)
@click.option('--ttmodel', type=str)
@click.option('--ttfeatures', type=str)
def cli(gene, ttmodel, ttfeatures):
    """
    Plots the blueprint
    """

    # Matrix for dotplot 
    fn = os.path.join(cancer_predictions_path, f'{gene}.model.{ttmodel}.features.{ttfeatures}.prediction.tsv.gz') 
    cancer_mat = pd.read_csv(fn, sep="\t")
    cancer_mat = cancer_mat[~(cancer_mat['aachange'].isna())].reset_index(drop=True)
    df_codon = load_saturation_cancer(cancer_mat, gene, shap_corrected=False)
    
    # Matrix for needleplot
    cancer_mat["Protein_position"] = cancer_mat.apply(lambda row: get_position(row), axis=1)
    
    # Observed mutations
    obs_mut.rename(columns={'mut': 'alt'}, inplace=True)
    obs_mut['chr'] = obs_mut['chr'].astype(str)
    obs_mut['pos'] = obs_mut['pos'].astype(int)
    df = create_observed_dataset(fn, gene, ttfeatures, obs_mut)
    
    # PLOT
    tracked_blueprint_all(gene, ttmodel, ttfeatures, df_codon, df, cancer_mat, show=True)
    plt.close()


if __name__ == '__main__':
    cli()