import click
import pandas as pd

pd.set_option('display.max_columns', 500)

from boostdm.globals import DRIVERS_PATH
from boostdm.oncotree import Oncotree
from boostdm.annotations.utils import encode_consequence_type, rectify_synonymous, rectify_missense, rectify_splicing
from boostdm.features import phylop, consequence_type, aachange, exon, ptms, clustl, hotmaps, smregions


drivers = pd.read_csv(DRIVERS_PATH, sep='\t')
dg = drivers.groupby('SYMBOL').agg({'CANCER_TYPE': list}).reset_index()
gene_ttype_map = dict(zip(dg['SYMBOL'].values, dg['CANCER_TYPE'].values))


def read_muts(path_data):

    muts = pd.read_csv(path_data, sep='\t')
    muts.rename(columns={'Chromosome': 'chr',
                         'Position': 'pos',
                         'Alternate': 'alt',
                         'Gene': 'ENSEMBL_GENE',
                         'Feature': 'ENSEMBL_TRANSCRIPT',
                         'Consequence': 'csqn_type',
                         'Symbol': 'gene',
                         'Canonical': 'CANONICAL'}, inplace=True)


    muts = muts[muts['CANONICAL'] == 'YES']
    if muts.shape[0] == 0:
        raise Exception('There are not mutations in the canonical transcript')

    muts['chr'] = muts['chr'].astype(str)
    muts['chr'] = muts['chr'].str.replace('chr', '')
    muts = muts[(muts['alt'].isin(['A', 'C', 'T', 'G']))]  # SNV

    return muts


def iter_tree(tree, ttype):
    """
    Given a query ttype and a tumor-type ontology (tree)
    generate all the tumor-types that are descendants of ttype in the tree.
    """

    res = []
    stack = [ttype]
    while len(stack) > 0:
        tt = stack.pop(0)
        res.append(tt)
        stack += tree.tree.get(tt, [])
    return res


def set_aachange(row):

    if row["Protein_position"] == "-":
        return "."
    else:
        pos = row["Protein_position"]
        wt_aa = row["Amino_acids"][0]
        mt_aa = row["Amino_acids"][-1]
        return str(wt_aa) + str(pos) + str(mt_aa)


def features(df, ttype, clustl_path, hotmaps_path, smregions_path):
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

    # tumor-type grouping:

    tree = Oncotree()
    related_ttypes = list(iter_tree(tree, ttype))

    # Add linear clusters
    clustl_global_data = pd.read_csv(clustl_path, sep='\t')
    clustl_cancer_data = clustl_global_data[clustl_global_data['CANCER_TYPE'].isin(related_ttypes)]
    clustl_cancer_data = clustl_cancer_data[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]
    clustl_global_data = clustl_global_data[['CHROMOSOME', '5_COORD', '3_COORD', 'SCORE']]
    df = clustl.add_feature(df, clustl_cancer_data, clustl_global_data)

    # Add 3D clusters
    hotmaps_global_data = pd.read_csv(hotmaps_path, sep='\t')
    hotmaps_cancer_data = hotmaps_global_data[hotmaps_global_data['CANCER_TYPE'].isin(related_ttypes)]
    hotmaps_cancer_data = hotmaps_cancer_data[['chromosome', 'pos']]
    hotmaps_global_data = hotmaps_global_data[['chromosome', 'pos']]
    df = hotmaps.add_feature(df, hotmaps_cancer_data, hotmaps_global_data)

    # run add_domains
    smregions_global_data = pd.read_csv(smregions_path, sep='\t')
    smregions_cancer_data= smregions_global_data[smregions_global_data['CANCER_TYPE'].isin(related_ttypes)]
    df = smregions.add_feature(df, smregions_cancer_data, smregions_global_data)

    return df


def build_table(mutations_file, tumor, path_clustl, path_hotmaps, path_smregions):

    # read mutations from the VEP output
    muts = read_muts(mutations_file)

    # Reset index
    muts.reset_index(inplace=True)

    # annotate mutations
    df = features(muts, tumor, path_clustl, path_hotmaps, path_smregions)

    # keep mutations with specified consequence type
    df = df[df['csqn_type'].isin(['synonymous', 'missense', 'nonsense', 'splicing'])]

    if df.shape[0] == 0:
        raise Exception("There are not 'synonymous', 'missense', 'nonsense', or 'splicing' mutations")

    df = encode_consequence_type(df)
    df = rectify_synonymous(df)
    df = rectify_missense(df)
    df = rectify_splicing(df)

    return df


@click.command()
@click.option('--gene', help="Gene symbol", type=str, required=True)
@click.option('--ttype', help="Tumor type acronym", type=str, required=True)
@click.option('--mutations', help="Path to the input mutations file from vep output",
              type=click.Path(exists=True), required=True)
@click.option('--clustl-group', help="CLUSTL clusters by ttype",
              type=click.Path(exists=True), required=True)
@click.option('--hotmaps-group', help="HotMAPs clusters by ttype",
              type=click.Path(exists=True),required=True)
@click.option('--smregions-group', help="SMRegions domains by ttype",
              type=click.Path(exists=True), required=True)
def cli(gene, ttype, mutations, clustl_group, hotmaps_group, smregions_group):
    df = build_table(mutations, ttype, clustl_group, hotmaps_group, smregions_group)
    df.to_csv(f'{gene}.{ttype}.annotated.tsv.gz', sep='\t', index=False, compression="gzip")

    # for testing only:
    # df = pd.DataFrame([])
    # df.to_csv(f'{gene}.{ttype}.annotated.tsv.gz', sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
