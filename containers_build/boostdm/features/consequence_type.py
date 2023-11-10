from functools import partial

from boostdm.globals import TABIX_FILE
from boostdm.vepreader import Tabix


# Consequence list taken from: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
CONSEQUENCES_LIST = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
]

CONSEQUENCES_DICT = {k: i for i, k in enumerate(CONSEQUENCES_LIST)}

AGGREGATION_DICT = {'synonymous_variant': 'synonymous',
                    'missense_variant': 'missense',
                    'stop_gained': 'nonsense',
                    'splice_donor_variant': 'splicing',
                    'splice_acceptor_variant': 'splicing',
                    'splice_region_variant': 'splicing'}


def get_csqn_type(chr_, pos, alt, gene, reader):

    for data in reader.get(chr_, pos, pos):
        alt_vep = (data[3] == alt)
        canonical_vep = (data[-4] == 'YES')  # impose canonical transcript
        correct_gene = (data[-7] == gene)  # skip few cases with antisense overlapping gene
        if alt_vep and canonical_vep and correct_gene:
            csqn = CONSEQUENCES_LIST[min([CONSEQUENCES_DICT[i] for i in data[7].split(',')])]
            return AGGREGATION_DICT.get(csqn, None)
    return None


def add_feature(df):

    with Tabix(TABIX_FILE) as reader:
        get_from_reader = partial(get_csqn_type, reader=reader)
        df['csqn_type'] = df.apply(lambda row: get_from_reader(str(row['chr']), int(row['pos']), row['alt'], row['gene']), axis=1)
    return df
