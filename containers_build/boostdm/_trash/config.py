"""Configuration object with all required constants of the project"""


import os


class PathConfig:

    intogen_data = os.environ['INTOGEN_DATASETS']
    drivers_path = os.environ['DRIVERS_PATH']
    cohorts_path = os.environ['COHORTS_PATH']
    cds = os.path.join(intogen_data, 'shared', 'cds.regions.gz')

    drivers = os.path.join(drivers_path)
    stats_cohort = os.path.join(cohorts_path)

    canonical_transcripts = os.path.join(intogen_data, 'shared', 'ensembl_canonical_transcripts.tsv')

    # Features
    phylop = os.environ['PHYLOP_FILE']
    pfam_domains = os.path.join(intogen_data, 'smregions', 'regions_pfam.tsv')
    ptms = os.path.join(intogen_data, 'ptms', 'info_functional_sites.json')


    # XGBoost params
    xgbparams = {
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

    run_folder = os.path.dirname(os.path.abspath(__file__))

    genes_json = os.path.join(intogen_data, 'shared', 'dictionary_datasets.json')

    columns_training = [
        'CLUSTL_SCORE', 'CLUSTL_cat_1', 'CLUSTL_cat_2',
        'HotMaps_cat_1', 'HotMaps_cat_2',
        'smRegions_cat_1', 'smRegions_cat_2',
        'PhyloP',
        'nmd',
        'Acetylation', 'Methylation', 'Phosphorylation', 'Regulatory_Site', 'Ubiquitination',
        'csqn_type_missense', 'csqn_type_nonsense', 'csqn_type_splicing', 'csqn_type_synonymous'
    ]

    columns_output = ['gene', 'ENSEMBL_TRANSCRIPT', 'ENSEMBL_GENE', 'chr', 'pos', 'alt', 'aachange'] + \
                     columns_training + \
                     ['selected_model_ttype', 'selected_model_gene', 'boostDM_score', 'boostDM_class']

    vep_tabix = os.path.join(intogen_data, 'shared', 'vep.tsv.bgz')

    splits = 50
    cv_fraction = 0.3

    min_test_size = 30
    auc_threshold = 0.8

    csqn_vals = ['missense', 'nonsense', 'splicing', 'synonymous']
    clustl_vals = [1, 2]
    hotmaps_vals = [1, 2]
    mnvs = os.path.join(intogen_data, 'boostDM', 'mnvs', 'mnvs.tsv.gz')

    # This option can be a number or None. In case of being None, the seed is not set.
    numpy_random_seed = 42
