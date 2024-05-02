import os

GENOME_BUILD = os.environ['GENOME_BUILD']
INTOGEN_DATASETS = os.environ['INTOGEN_DATASETS']
BOOSTDM_DATASETS = os.environ['BOOSTDM_DATASETS']

COHORTS_PATH = os.path.join(INTOGEN_DATASETS, 'cohorts.tsv')
DRIVERS_PATH = os.path.join(INTOGEN_DATASETS, 'drivers.tsv')
PHYLOP_FILE = os.path.join(BOOSTDM_DATASETS, 'hg38.phyloP100way.bw')
MNVS_FILE = os.path.join(INTOGEN_DATASETS, 'steps', 'boostDM', 'mnvs.tsv.gz')
MANE_TRANSCRIPTS_FILE = os.path.join(BOOSTDM_DATASETS, 'saturation', 'cds-25spli.regions.gz')
TABIX_FILE = os.path.join(BOOSTDM_DATASETS, 'shared', 'vep.tsv.gz')
PTMS_FILE = os.path.join(BOOSTDM_DATASETS, 'ptms', 'info_functional_sites.json')
PFAM_DOMAINS_FILE = os.path.join(BOOSTDM_DATASETS, 'regions_pfam.tsv')
ONCOTREE_PATH = os.path.join(BOOSTDM_DATASETS, 'shared', 'tree.tsv')

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
    'splice_donor_5th_base_variant',
    'splice_region_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
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
    'intergenic_variant',
    'sequence_variant'
]

CONSEQUENCES_DICT = {k: i for i, k in enumerate(CONSEQUENCES_LIST)}

AGGREGATION_DICT = {'synonymous_variant': 'synonymous',
                    'missense_variant': 'missense',
                    'stop_gained': 'nonsense',
                    'stop_lost': 'nonsense',
                    'start_lost': 'nonsense',
                    'splice_donor_variant': 'splicing',
                    'splice_acceptor_variant': 'splicing',
                    'splice_region_variant': 'splicing',
                    'intron_variant': 'splicing'}

COLUMNS_TRAINING = [
        'CLUSTL', 'HotMaps', 'smRegions', 'PhyloP',
        'nmd', 
        'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination',
        'csqn_type_missense', 'csqn_type_nonsense', 'csqn_type_splicing', 'csqn_type_synonymous'
    ]

COLUMNS_OUTPUT = ['gene', 'ENSEMBL_TRANSCRIPT', 'ENSEMBL_GENE', 'chr', 'pos', 'alt', 'aachange'] + \
                 COLUMNS_TRAINING + \
                 ['selected_model_ttype'] + \
                 ['boostDM_score', 'boostDM_class']

COLUMNS_SHAP = [f'shap_{x}' for x in COLUMNS_TRAINING]


# XGBoost params
XGB_PARAMS = {
        "objective": "binary:logistic",
        "reg_lambda": 1,
        "random_state": 42,
        "scale_pos_weight": 1,
        "subsample": 1.0,
        "reg_alpha": 0,
        "max_delta_step": 0,
        "min_child_weight": 1,
        "learning_rate": 1e-03,
        "colsample_bylevel": 1.0,
        "gamma": 0,
        "colsample_bytree": 1.0,
        "booster": "gbtree",
        "max_depth": 4,
        "silent": True,
        "seed": 21
}


# run: output_20230323

# DISCOVERY_TIERS = (0, 0.2, 0.4, 0.6, 0.8)
# MUTATION_TIERS = (40, 30, 20, 10, 5)
# FSCORE_THRESHOLD = 0.8


# run: output_20230503

DISCOVERY_TIERS = (0, 0.5, 0.75)
MUTATION_TIERS = (50, 30, 0)
FSCORE_THRESHOLD = 0.8


# consensus merge parameter

SYSTEMATIC_BIAS = 2.3
