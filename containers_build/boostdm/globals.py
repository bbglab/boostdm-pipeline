import os

INTOGEN_DATASETS = os.environ['INTOGEN_DATASETS']

COHORTS_PATH = os.path.join(INTOGEN_DATASETS, 'cohorts.tsv')
DRIVERS_PATH = os.path.join(INTOGEN_DATASETS, 'drivers.tsv')

ONCOTREE_PATH = os.path.join(INTOGEN_DATASETS, "debug", "oncotree", "tree_cancer_types.json")

CDS_FILE = os.path.join(INTOGEN_DATASETS, 'debug', 'shared', 'cds.regions.gz')
PHYLOP_FILE = os.path.join(INTOGEN_DATASETS, 'debug', 'phylop', 'hg38.phyloP100way.bw')
TABIX_FILE = os.path.join(INTOGEN_DATASETS, 'debug', 'shared', 'vep.tsv.bgz')
PTMS_FILE = os.path.join(INTOGEN_DATASETS, 'debug', 'ptms', 'info_functional_sites.json')
PFAM_DOMAINS_FILE = os.path.join(INTOGEN_DATASETS, 'debug', 'smregions', 'regions_pfam.tsv')


COLUMNS_TRAINING = [
        'CLUSTL_SCORE', 'CLUSTL_cat_1',
        'HotMaps_cat_1',
        'smRegions_cat_1',
        'PhyloP',
        'nmd',
        'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination',
        'csqn_type_missense', 'csqn_type_nonsense', 'csqn_type_splicing', 'csqn_type_synonymous'
    ]  
#  "cat_2" features for CLUSTL, Hotmaps and smRegions have been removed


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
        "subsample": 0.7,
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


# used in run_20230725

# DISCOVERY_TIERS = (0, 0.2, 0.4, 0.6, 0.8)
# MUTATION_TIERS = (40, 30, 20, 10, 5)
# FSCORE_THRESHOLD = 0.8


# used in run_20230802

DISCOVERY_TIERS = (0, 0.5, 0.75)
MUTATION_TIERS = (50, 30, 0)
FSCORE_THRESHOLD = 0.8


# consensus merge parameter

SYSTEMATIC_BIAS = 2.3
