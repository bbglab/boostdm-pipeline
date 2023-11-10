import os
import gzip
import pickle
import click
import pandas as pd

from boostdm.config import PathConfig
from boostdm.oncotree import Oncotree

from boostdm.annotations.utils import encoding, rectify_synonymous, rectify_missense, rectify_splicing
from boostdm.features import phylop, consequence_type, aachange, exon, ptms, clustl, hotmaps, smregions


# TODO do not use VEP as input, read the information from the VEP file if not provided
def read_muts(path_data):

    muts = pd.read_csv(path_data, sep="\t")

    muts = muts[muts['Canonical'] == 'YES']
    if muts.shape[0] == 0:
        raise Exception('There are not mutations in canonical transcripts')

    muts.rename(columns={
        'Chromosome': 'chr',
        'Position': 'pos',
        'Alternate': 'alt',
        'Symbol': 'gene',
        'Consequence': 'csqn_type',
        'Gene': 'ENSEMBL_GENE',
        'Feature': 'ENSEMBL_TRANSCRIPT'
    }, inplace=True)

    muts = muts[(muts['alt'].isin(['A', 'C', 'T', 'G']))]  # SNV
    muts['chr'] = muts['chr'].astype(str)
    muts['chr'] = muts['chr'].str.replace('chr', '')

    return muts


def iter_tree(tree, ttype):
    yield ttype
    if ttype in tree.tree:
        for child in tree.tree[ttype]:
            yield iter_tree(tree, child)


def filter_drivers(muts):

    # Filter out mutations in genes not annotated as drivers
    df_drivers = pd.read_csv(PathConfig.drivers, sep="\t")
    drivers = df_drivers["SYMBOL"].unique()
    return muts[muts["gene"].isin(drivers)]


def set_aachange(row):

    if row["Protein_position"] == "-":
        return "."
    else:
        pos = row["Protein_position"]
        wt_aa = row["Amino_acids"][0]
        mt_aa = row["Amino_acids"][-1]
        return str(wt_aa) + str(pos) + str(mt_aa)


def features(df, tumor, drivers_summary, clustl_path, hotmaps_path, smregions_path):
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

    tree = Oncotree()
    related_ttypes = [ttype for ttype in iter_tree(tree, tumor)]

    # Add linear clusters
    clustl_global_data = pd.read_csv(clustl_path, sep='\t')
    clustl_cancer_data = clustl_global_data[clustl_global_data['CANCER_TYPE'].isin(related_ttypes)]

    # TODO .drop_duplicates() needed ?
    clustl_cancer_data = clustl_cancer_data[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]
    # TODO .drop_duplicates() needed ?
    clustl_global_data = clustl_global_data[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]
    df = clustl.add_feature(df, clustl_cancer_data, clustl_global_data)

    # Add 3D clusters
    hotmaps_global_data = pd.read_csv(hotmaps_path, sep='\t')
    hotmaps_cancer_data = hotmaps_global_data[hotmaps_global_data['CANCER_TYPE'].isin(related_ttypes)]

    # TODO .drop_duplicates() needed
    hotmaps_cancer_data = hotmaps_cancer_data[['chromosome', 'pos']]
    # TODO .drop_duplicates() needed ?
    hotmaps_global_data = hotmaps_global_data[['chromosome', 'pos']]
    df = hotmaps.add_feature(df, hotmaps_cancer_data, hotmaps_global_data)

    # add role
    df_role = pd.read_csv(drivers_summary, sep='\t')
    df_role.rename(columns={'SYMBOL': 'gene',
                            'ROLE': 'role'},
                   inplace=True)
    df = df.merge(df_role[['gene', 'role']].drop_duplicates())

    # run add_domains
    smregions_global_data = pd.read_csv(smregions_path, sep='\t')
    smregions_cancer_data= smregions_global_data[smregions_global_data['CANCER_TYPE'].isin(related_ttypes)]
    df = smregions.add_feature(df, smregions_cancer_data, smregions_global_data)

    return df


def build_table(mutations_file, ttype, drivers_summary_path, path_clustl, path_hotmaps, path_smregions):

    # read mutations from the VEP output
    muts = read_muts(mutations_file)

    # Reset index
    muts.reset_index(drop=True, inplace=True)

    # climb the ontology until the top level containing ttype
    # the resulting tt will be the level used for feature annotation
    tree = Oncotree()
    tt = ttype
    parent_ttype = tree.fetch_parent_ttype(tt)
    while parent_ttype not in ['STOP', None]:
        tt = parent_ttype
        parent_ttype = tree.fetch_parent_ttype(tt)

    # filter mutations in drivers
    muts_drivers = filter_drivers(muts)
    if muts_drivers.shape[0] == 0:
        raise Exception("There are not driver mutations")

    # annotate mutations
    mut_features = features(muts_drivers, tt, drivers_summary_path, path_clustl, path_hotmaps, path_smregions)

    # keep mutations with specified consequence type
    mut_features = mut_features[mut_features['csqn_type'].isin(['synonymous',
                                                                'missense',
                                                                'nonsense',
                                                                'splicing'])]

    if mut_features.shape[0] == 0:
        raise Exception("There are not 'synonymous', 'missense', 'nonsense', or 'splicing' mutations")

    # encode the mutations
    df = encoding(mut_features)

    # run without pancancer feature annotations
    df = df.drop(columns=['CLUSTL_cat_2',
                          'HotMaps_cat_2',
                          'smRegions_cat_2'])

    # rectify
    df = rectify_synonymous(df)
    df = rectify_missense(df)
    df = rectify_splicing(df)

    # feature origin: ttype-type that has been used to derive the features
    df['feature_origin'] = tt

    return df


@click.command()
@click.option(
    '--mutations',
    help="Path to the input mutations file from vep output",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--drivers-summary',
    help="Drivers summary file",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--clustl-group',
    help="ClustL clusters by ttype",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--hotmaps-group',
    help="HotMAPs clusters by ttype",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--smregions-group',
    help="SMRegions domains by ttype",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--ttype',
    help="Name of the ttype to be used for the annotation of the features",
    type=str,
    default=None
)
@click.option(
    '--output',
    help="Output file to write the annotated mutations",
    type=click.Path()
)
def cli(mutations, ttype, drivers_summary, clustl_group, hotmaps_group, smregions_group, output):

    df = build_table(mutations, ttype, drivers_summary, clustl_group, hotmaps_group, smregions_group)

    # save them
    df.to_csv(output, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
