
import click
import pandas as pd

from boostdm.globals import PFAM_DOMAINS_FILE


def add_feature(df, specific_df, global_df):

    # TODO: implement as an apply
    
    isinmotifs = []
    significant_motif = []

    # classify mutations
    # TODO get rid of iterrows
    for i, row in df.iterrows():
        mut_in_motif = 3  # it was 0, changed to 3 for consistence with other methods
        motif_sig = ''
        symbol = row['gene']
        
        # select only the gene symbol to speed things up
        motifs_in_trans = specific_df[specific_df['SYMBOL'] == symbol]
        if len(motifs_in_trans) > 0:
            for ix, line in motifs_in_trans.iterrows():
                if (int(line['START']) < row['pos']) and (row['pos'] < int(line['STOP'])) and (row['csqn_type'] != "synonymous_variant"):
                    mut_in_motif = 1
                    motif_sig = '{};{}'.format(line['ELEMENT2'], motif_sig)

        if mut_in_motif == 0:
            motifs_in_trans = global_df[global_df['SYMBOL'] == symbol]
            if len(motifs_in_trans) > 0:
                for ix, line in motifs_in_trans.iterrows():
                    if (int(line['START']) < row['pos']) and (row['pos'] < int(line['STOP'])) and (
                            row['csqn_type'] != "synonymous_variant"):
                        mut_in_motif = 2
                        motif_sig = '{};{}'.format(line['ELEMENT2'], motif_sig)

        isinmotifs.append(mut_in_motif)
        significant_motif.append(motif_sig)

    df['motif'] = significant_motif
    df['smRegions_cat'] = isinmotifs
    return df


def motifs(significant_regions):

    pfam_regions = pd.read_csv(PFAM_DOMAINS_FILE,
                               sep='\t', low_memory=False,
                               names=['CHROMOSOME', 'START', 'STOP', 'STRAND', 'ELEMENT', 'ELEMENT2', 'SYMBOL'])

    # select only significant regions!
    selected_trans_pfam = pfam_regions[pfam_regions['ELEMENT'].isin(significant_regions)]

    return selected_trans_pfam


def generate(files_list, qval_thresh=0.1):
    """group a set of output files into a single dataframe"""

    df = []
    for input_file in files_list:
        df_dom = pd.read_csv(input_file, sep='\t')
        df.append(df_dom[(df_dom['Q_VALUE'] < qval_thresh) & ((df_dom['OBSERVED_REGION'] / df_dom['MEAN_SIMULATED']) > 1)])
    df = pd.concat(df, axis=0)

    significant_regions = df['REGION'].unique().tolist()
    if len(significant_regions) > 0:
        df = motifs(significant_regions)
    else:
        df = pd.DataFrame([],columns=['CHROMOSOME', 'START', 'STOP', 'STRAND', 'ELEMENT', 'ELEMENT2', 'SYMBOL'])

    return df


@click.command()
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--threshold', default=0.1, type=float, help='Qvalue threshold')
@click.argument('files', nargs=-1)
def cli(files, output, threshold):

    df = generate(files, qval_thresh=threshold)

    df.to_csv(output, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    cli()
