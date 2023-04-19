# Imports
# -------
import itertools
import os
import warnings

import click
import numpy as np
import pandas as pd
import pybedtools
import tqdm

from boostdm import BoostDMError
from boostdm.annotations.utils import encoding, rectify_synonymous, rectify_missense, rectify_splicing
from boostdm.globals import CANONICAL_TRANSCRIPTS_FILE, MNVS_FILE
from boostdm.mutrate_reader import MutrateReader
from boostdm.oncotree import Oncotree
from boostdm.features import phylop, consequence_type, aachange, exon, ptms, clustl, hotmaps, smregions, dndscv
from boostdm.passengers import retrieve_exons, randomize

COLUMNS = ['sampleID', 'chr', 'start', 'end', 'ref', 'alt', 'gene']

CB = dict(zip(list('ACGT'), list('TGCA')))


# Utils
# -----

def mut_key_generator():

    subs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
    for s in sorted(subs):
        for c in sorted(itertools.product({'A', 'C', 'G', 'T'}, repeat=2)):
            yield tuple([s, ''.join(c)])


def purine_to_pyrimidine(m):
    if m[0][0] not in list('CT'):
        t1 = CB[m[0][0]] + CB[m[0][1]]
        t2 = CB[m[1][1]] + CB[m[1][0]]
        return tuple([t1, t2])
    else:
        return m


mut_keys = [k for k in mut_key_generator()]


def retrieve_expectation(exp_dict, v):
    sample = v['sampleID']
    ref3_cod = v['ref3_cod']
    mut3_cod = v['mut3_cod']
    key = tuple([ref3_cod[1] + mut3_cod[1], ref3_cod[0] + ref3_cod[2]])
    key = purine_to_pyrimidine(key)
    if sample in exp_dict:
        return exp_dict[sample][mut_keys.index(key)]
    else:
        return 0.0


def set_string_chr(row):
    try:
        return str(int(row["chr"]))
    except:
        return str(row["chr"])


# Pipeline
# --------


def oncotree_sisters(cohort):
    """Generator of cohorts belonging to the same tumor type as 'cohort'"""
    # TODO re-implement function in terms of oncotree
    tree = Oncotree()
    parent = tree.fetch_parent_cohort(cohort)
    cohorts = tree.get_cohorts(parent)
    return cohorts


def features(df, drivers_summary,
             clustl_path, clustl_group_path,
             hotmaps_path, hotmaps_group_path,
             smregions_path, smregions_group_path):
    """add complete set of features"""

    # add PhyloP score
    df = phylop.add_feature(df)

    # add consequence type
    df = consequence_type.add_feature(df)  # from vep

    # add amino-acid change
    df = aachange.add_feature(df)  # from vep

    # add NMD annotation:
    df = exon.add_feature(df)  # from vep

    # add PTMs:
    df = ptms.add_ptms(df)  # from phosphosite-plus

    # Add linear clusters
    clustl_cohort_data = clustl.generate([clustl_path])
    clustl_global_data = pd.read_csv(clustl_group_path, sep='\t')
    clustl_global_data = clustl_global_data[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]  # TODO .drop_duplicates() needed ?
    df = clustl.add_feature(df, clustl_cohort_data, clustl_global_data)

    # Add 3D clusters
    hotmaps_cohort_data = hotmaps.generate([hotmaps_path])
    hotmaps_global_data = pd.read_csv(hotmaps_group_path, sep='\t')
    hotmaps_global_data = hotmaps_global_data[['chromosome', 'pos']]   # TODO .drop_duplicates() needed ?
    df = hotmaps.add_feature(df, hotmaps_cohort_data, hotmaps_global_data)

    # add role
    df_role = pd.read_csv(drivers_summary, sep='\t')
    df_role.rename(columns={'SYMBOL': 'gene', 'ROLE': 'role'}, inplace=True)
    df = df.merge(df_role[['gene', 'role']].drop_duplicates())

    # run add_domains
    smregions_cohort_data = smregions.generate([smregions_path])
    smregions_global_data = pd.read_csv(smregions_group_path, sep='\t')
    smregions_global_data = smregions_global_data  # TODO .drop_duplicates() needed ?
    df = smregions.add_feature(df, smregions_cohort_data, smregions_global_data)

    return df


def mnvs_to_remove():
    """
    Produces a set of SNVs that must be removed,
    just because they may be MNVs in disguise
    """

    mnvs_df = pd.read_csv(MNVS_FILE, sep='\t')
    mnvs_df = mnvs_df[mnvs_df['MNV']]
    mnvs_list = list(zip(mnvs_df.SYMBOL, mnvs_df.sample_id, mnvs_df.COHORT, mnvs_df.pos.str.split(',')))
    s = set([(x[0], x[1], x[2], int(a)) for x in mnvs_list for a in x[3]])
    return s


# Functions
# ---------


def retrieve_cds():

    """Returns dataframe with Canonical Transcript regions"""

    cds_df = pd.read_csv(CANONICAL_TRANSCRIPTS_FILE, sep='\t', header=None, compression='gzip', low_memory=False, skiprows=1)
    cds_df = cds_df[[0, 1, 2, 6]].copy()
    cds_df.columns = ['chr', 'start', 'end', 'gene']
    return cds_df


def intersect_region_mutations(cds, pos):

    """
    Returns:
         data frame with true positives, including sample and gene annotation
         dictionary of counts per sample-gene
    """

    pos = pos[['chr', 'start', 'end', 'ref', 'alt', 'gene']]
    cds_bed = pybedtools.BedTool.from_dataframe(cds)
    pos_bed = pybedtools.BedTool.from_dataframe(pos)
    drivers = pos_bed.intersect(cds_bed, wao=True).to_dataframe(
        usecols=[0, 2, 3, 4, 5],
        names=['chr', 'pos', 'ref', 'alt', 'gene']
    )

    return drivers


def load_drivers(cohort, df):

    # ttype = df[df['COHORT'] == cohort]['CANCER_TYPE'].unique()[0]
    # TODO replace this and get a list of useful cancer types:
    # beware that some cohorts might not have drivers, so the step above
    # will fail but the sister cohorts still can have drivers, so we
    # must use a different mapping
    cohort_list = oncotree_sisters(cohort)
    df_drivers_summary = df[df['COHORT'].isin(cohort_list)]

    drivers = df_drivers_summary['SYMBOL'].unique()

    return drivers


def load_mutations(dndscv_annotated, drivers, cohort):
    df = pd.read_csv(dndscv_annotated, sep='\t')  # os.path.join(dndscv_path, f'{cohort}.annotmuts.gz')

    df = df[df['gene'].isin(drivers)]
    if len(df) == 0:
        raise BoostDMError('Run failed: empty dataframe')

    # TODO how to compute the MNVS?
    # Filter out possible double substitutions
    s = mnvs_to_remove()

    # TODO cohort param is only required for the mnvs removal step
    df['cohort'] = cohort
    df['tuple'] = df.apply(lambda x: (x['gene'], x['sampleID'], x['cohort'], x['pos']), axis=1)
    df = df[~df['tuple'].isin(s)]
    del df['tuple']

    return df


def compute_expected(df, drivers, dict_genes):

    list_dfs = []
    for gene in tqdm.tqdm(drivers):

        gene_dict_total = None

        try:
            gene_dict_total = dict_genes[gene]
        except KeyError:
            warnings.warn(f'mutrate is not available for gene={gene}', UserWarning)
            continue

        dg = df[(df['gene'] == gene) & (df['impact'] != "no-SNV")]
        if dg.shape[0] == 0:
            continue
        dg['expect'] = dg.apply(lambda v: retrieve_expectation(gene_dict_total, v), axis=1)
        list_dfs.append(dg)

    return pd.concat(list_dfs)


def build_positive_set(df_expect):
    cds = retrieve_cds()
    pos = intersect_region_mutations(cds, df_expect)
    pos['response'] = 1
    return pos


def build_negative_set(counts_per_gene, dict_genes, n_splits=50):

    neg_dict = {'chr': [], 'pos': [], 'ref': [], 'alt': [], 'gene': []}

    for gene in tqdm.tqdm(counts_per_gene):

        # get mutrates of the entire cohort
        mutrates = dict_genes[gene]
        if mutrates is None:
            continue
        n_channels = len(list(mutrates.values())[0])
        total_mutrate = np.zeros(n_channels)
        for sample in mutrates:
            total_mutrate += np.array(mutrates[sample])

        # get exon information per gene
        chrom, cds, segments = retrieve_exons(gene)
        for chr_, mut in randomize(total_mutrate, chrom, cds, segments, n_splits * counts_per_gene[gene]):
            neg_dict['chr'].append(chr_)
            neg_dict['pos'].append(mut.pos)
            neg_dict['ref'].append(mut.ref_triplet[1])
            neg_dict['alt'].append(mut.alt)
            neg_dict['gene'].append(gene)

    return pd.DataFrame(neg_dict)


def add_passenger_mutations(df, dict_genes, n_splits):

    # create dict of counts per gene
    g = df.groupby(['gene']).size()
    gene_count = g.to_dict()

    neg = build_negative_set(gene_count, dict_genes, n_splits=n_splits)
    neg['response'] = 0
    return neg


def build_table(cohort, drivers_summary_file, dndscv_file, dndscv_annotated_file,
                mutrate, clustl_file, clustl_group_file,
                hotmaps_file, hotmaps_group_file,
                smregions_file, smregions_group_file, xs_thresh=0.85, n_splits=50):

    df_drivers_summary = pd.read_csv(drivers_summary_file, sep='\t')

    drivers = load_drivers(cohort, df_drivers_summary)
    if len(drivers) == 0:
        raise BoostDMError('Run failed: no drivers for this cohort')

    df = load_mutations(dndscv_annotated_file, drivers, cohort)

    dict_genes = {}

    # Read mutrate from the tarfile
    # FIXME load as a single JSON when changes are implemented in intogen+
    if os.path.exists(mutrate):
        mr = MutrateReader(mutrate)
        for gene in tqdm.tqdm(drivers):
            dict_genes[gene] = mr.load(gene)

    df_expect = compute_expected(df, drivers, dict_genes)

    df_expect = dndscv.filter(df_expect, dndscv_file, xs_thresh=xs_thresh)

    # define columns of the output dataset
    df_expect.rename(columns={'mut': 'alt'}, inplace=True)
    df_expect.rename(columns={'pos': 'end'}, inplace=True)
    df_expect['start'] = df_expect['end'].apply(lambda x: int(x) - 1)
    df_expect['end'] = df_expect['end'].apply(lambda x: int(x))

    df_expect = df_expect[COLUMNS]

    pos = build_positive_set(df_expect)
    neg = add_passenger_mutations(pos, dict_genes, n_splits)

    df = pd.concat([pos, neg], sort=True)

    # Format columns
    df['chr'] = df.apply(lambda row: set_string_chr(row), axis=1)
    df['pos'] = df.apply(lambda row: int(row['pos']), axis=1)

    # Add features
    df = features(df, drivers_summary_file,
             clustl_file, clustl_group_file,
             hotmaps_file, hotmaps_group_file,
             smregions_file, smregions_group_file)
    df = encoding(df)
    df = df[(df['csqn_type_synonymous'] != 1) | (df['response'] != 1)]
    df = rectify_synonymous(df)
    df = rectify_missense(df)
    df = rectify_splicing(df)

    return df


@click.command()
@click.option('--cohort', type=str)
@click.option('--drivers-summary', 'summary', type=click.Path(exists=True), help='Drivers summary from IntOGen')
@click.option('--dndscv-path', 'dndscv_path', type=click.Path(exists=True), help='Cohort dNdsCV out')
@click.option('--dndscv-annotmuts-path', 'dnds_muts_path', type=click.Path(exists=True), help='Cohort dNdsCV annotmuts out')
@click.option('--mutrate-path', 'mutrate_path', type=click.Path(exists=True), help='Cohort mutrate out')
@click.option('--clustl-path', 'clustl_path', type=click.Path(exists=True), help='Cohort OncodriveCLUSTL out')
@click.option('--clustl-group-path', 'clustl_group_path', type=click.Path(exists=True), help='Combined OncodriveCLUSTL out')
@click.option('--hotmaps-path', 'hotmaps_path', type=click.Path(exists=True), help='Cohort HotMAPS out')
@click.option('--hotmaps-group-path', 'hotmaps_group_path', type=click.Path(exists=True), help='Combined HotMAPS out')
@click.option('--smregions-path', 'smregions_path', type=click.Path(exists=True), help='Cohort smregions out')
@click.option('--smregions-group-path', 'smregions_group_path', type=click.Path(exists=True), help='Combined smregions out')
@click.option('--out', type=click.Path())
@click.option('--seed', type=int, default=None)
@click.option('--splits', type=int, default=50)
@click.option('--threshold', type=float, default=0.85)
def cli(cohort, summary, dndscv_path, dnds_muts_path, mutrate_path, clustl_path, clustl_group_path,
        hotmaps_path, hotmaps_group_path, smregions_path, smregions_group_path, out, seed, splits, threshold):
    """build raw mutations table"""

    np.random.seed(seed)

    df = build_table(cohort, summary, dndscv_path, dnds_muts_path, mutrate_path, clustl_path,
                     clustl_group_path, hotmaps_path, hotmaps_group_path, smregions_path, smregions_group_path,
                     n_splits=splits, xs_thresh=threshold)

    df.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    cli()
