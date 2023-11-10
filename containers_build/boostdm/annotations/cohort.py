# Imports
# -------
import itertools
import click
import json
import numpy as np
import pandas as pd
import pybedtools
import tqdm

from boostdm import BoostDMError
from boostdm.annotations.utils import encode_consequence_type, rectify_synonymous, rectify_missense, rectify_splicing
from boostdm.globals import CANONICAL_TRANSCRIPTS_FILE, MNVS_FILE, COHORTS_PATH, DRIVERS_PATH
from boostdm.oncotree import Oncotree
from boostdm.features import phylop, consequence_type, aachange, exon, ptms, clustl, hotmaps, smregions, dndscv
from boostdm.passengers import retrieve_exons, randomize

COLUMNS = ['sampleID', 'chr', 'start', 'end', 'ref', 'alt', 'gene']

CB = dict(zip(list('ACGT'), list('TGCA')))


def load_ttypes_map(file):

    df = pd.read_csv(file, sep='\t')
    df = df[['COHORT', 'CANCER_TYPE']]
    return df.set_index('COHORT').to_dict()['CANCER_TYPE']


ttype_map = load_ttypes_map(COHORTS_PATH)


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
    
    """Generator of cohorts belonging to the same ttype type as 'cohort'"""
    # TODO re-implement function in terms of oncotree
    tree = Oncotree()
    parent = tree.fetch_parent_cohort(cohort)
    cohorts = tree.get_cohorts(parent)
    return cohorts


def features(df, cohort, clustl_group_path, hotmaps_group_path, smregions_group_path):

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

    # add tumor-type level and pan-cancer level cluster features
    ttype = ttype_map[cohort]

    # Add linear clusters

    clustl_global_data = pd.read_csv(clustl_group_path, sep='\t')
    clustl_ttype_data = clustl_global_data[clustl_global_data['CANCER_TYPE'] == ttype]
    # clustl_global_data = clustl_global_data[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]
    df = clustl.add_feature(df, clustl_ttype_data, clustl_global_data)

    # Add 3D clusters
    hotmaps_global_data = pd.read_csv(hotmaps_group_path, sep='\t')
    hotmaps_ttype_data = hotmaps_global_data[hotmaps_global_data['CANCER_TYPE'] == ttype]
    # hotmaps_global_data = hotmaps_global_data[['chromosome', 'pos']]
    df = hotmaps.add_feature(df, hotmaps_ttype_data, hotmaps_global_data)

    # add role
    # df_role = pd.read_csv(DRIVERS_PATH, sep='\t')
    # df_role.rename(columns={'SYMBOL': 'gene', 'ROLE': 'role'}, inplace=True)
    # df = df.merge(df_role[['gene', 'role']].drop_duplicates())

    # run add_domains
    smregions_global_data = pd.read_csv(smregions_group_path, sep='\t')
    smregions_ttype_data = smregions_global_data[smregions_global_data['CANCER_TYPE'] == ttype]
    df = smregions.add_feature(df, smregions_ttype_data, smregions_global_data)

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


def retrieve_transcript():

    """Returns dataframe with canonical transcript regions"""

    canonical_transcript_df = pd.read_csv(CANONICAL_TRANSCRIPTS_FILE,
                         sep='\t', header=None, compression='gzip', low_memory=False, skiprows=1)

    # TODO: verify the columns we are selecting are the right ones
    
    canonical_transcript_df = canonical_transcript_df[[0, 1, 2, 6]].copy()
    canonical_transcript_df.columns = ['chr', 'start', 'end', 'gene']
    return canonical_transcript_df


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


def initialize_trainset(df, drivers):
    list_dfs = []
    for gene in tqdm.tqdm(drivers):
        dg = df[(df['gene'] == gene) & (df['impact'] != "no-SNV")]
        if dg.shape[0] == 0:
            continue
        list_dfs.append(dg)
    return pd.concat(list_dfs)


def build_positive_set(df_expect):

    canonical_transcript = retrieve_transcript()
    pos = intersect_region_mutations(canonical_transcript, df_expect)
    pos['response'] = 1
    return pos


def build_negative_set(counts_per_gene, mutrate, n_splits=50):

    neg_dict = {'chr': [], 'pos': [], 'ref': [], 'alt': [], 'gene': []}

    for gene in tqdm.tqdm(counts_per_gene):

        if len(mutrate) == 0:
            continue

        # get exon information per gene
        chrom, cds, segments = retrieve_exons(gene)
        for chr_, mut in randomize(mutrate, chrom, cds, segments, n_splits * counts_per_gene[gene]):
            neg_dict['chr'].append(chr_)
            neg_dict['pos'].append(mut.pos)
            neg_dict['ref'].append(mut.ref_triplet[1])
            neg_dict['alt'].append(mut.alt)
            neg_dict['gene'].append(gene)

    return pd.DataFrame(neg_dict)


def add_passenger_mutations(df, mutrate, n_splits):

    # create dict of counts per gene
    g = df.groupby(['gene']).size()
    gene_count = g.to_dict()

    neg = build_negative_set(gene_count, mutrate, n_splits=n_splits)
    neg['response'] = 0
    return neg


def build_table(cohort, dndscv_file, dndscv_annotated_file,
                mutrate_fn, clustl_group_file, hotmaps_group_file, smregions_group_file,
                xs_thresh=0.85, n_splits=50):

    df_drivers_summary = pd.read_csv(DRIVERS_PATH, sep='\t')

    drivers = load_drivers(cohort, df_drivers_summary)
    if len(drivers) == 0:
        raise BoostDMError('Run failed: no drivers for this cohort')

    df = load_mutations(dndscv_annotated_file, drivers, cohort)

    # Read mutrate

    with open(mutrate_fn, 'rt') as f:
        mutrate = json.load(f)

    initial = initialize_trainset(df, drivers)
    initial = dndscv.filter(initial, dndscv_file, xs_thresh=xs_thresh)
    initial.rename(columns={'mut': 'alt'}, inplace=True)
    initial.rename(columns={'pos': 'end'}, inplace=True)
    initial['start'] = initial['end'].apply(lambda x: int(x) - 1)
    initial['end'] = initial['end'].apply(lambda x: int(x))
    initial = initial[COLUMNS]

    # positive set
    pos = build_positive_set(initial)

    # negative set

    neg = add_passenger_mutations(pos, mutrate, n_splits)

    # join pos + neg set

    df = pd.concat([pos, neg], sort=True)

    # Format columns
    df['chr'] = df.apply(lambda row: set_string_chr(row), axis=1)
    df['pos'] = df.apply(lambda row: int(row['pos']), axis=1)

    # Reset index
    df.reset_index(drop=True, inplace=True)

    # Add features
    df = features(df, cohort, clustl_group_file, hotmaps_group_file, smregions_group_file)
    df = encode_consequence_type(df)
    df = df[(df['csqn_type_synonymous'] != 1) | (df['response'] != 1)]
    df = rectify_synonymous(df)
    df = rectify_missense(df)
    df = rectify_splicing(df)
    
    return df


@click.command()
@click.option('--cohort', type=str)
@click.option('--dndscv-path', 'dndscv_path', type=click.Path(exists=True), help='Cohort dNdsCV out')
@click.option('--dndscv-annotmuts-path', 'dnds_muts_path', type=click.Path(exists=True), help='Cohort dNdsCV annotmuts out')
@click.option('--mutrate-path', 'mutrate_path', type=click.Path(exists=True), help='Cohort mutrate out')
@click.option('--clustl-group-path', 'clustl_group_path', type=click.Path(exists=True), help='Combined OncodriveCLUSTL out')
@click.option('--hotmaps-group-path', 'hotmaps_group_path', type=click.Path(exists=True), help='Combined HotMAPS out')
@click.option('--smregions-group-path', 'smregions_group_path', type=click.Path(exists=True), help='Combined smregions out')
@click.option('--out', type=click.Path())
@click.option('--seed', type=int, default=None)
@click.option('--splits', type=int, default=50)
@click.option('--threshold', type=float, default=0.85)
def cli(cohort, dndscv_path, dnds_muts_path, mutrate_path, clustl_group_path,
        hotmaps_group_path, smregions_group_path, out, seed, splits, threshold):
    """build raw mutations table"""

    np.random.seed(seed)

    df = build_table(cohort, dndscv_path, dnds_muts_path, mutrate_path,
                     clustl_group_path, hotmaps_group_path, smregions_group_path,
                     n_splits=splits, xs_thresh=threshold)

    df.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    cli()
