#!/usr/bin/env nextflow

/* COHORTS SUMMART */

COHORTS_SUMMARY = Channel.fromPath("${INTOGEN_DATASETS}/cohorts.tsv")
COHORTS_SUMMARY.into{ COHORTS_SUMMARY1; COHORTS_SUMMARY2; COHORTS_SUMMARY3;  COHORTS_SUMMARY4}

/* DRIVERS SUMMARY */

DRIVERS_SUMMARY = Channel.fromPath("${INTOGEN_DATASETS}/drivers.tsv")
DRIVERS_SUMMARY_V = DRIVERS_SUMMARY.first()

/* SATURATION ANALYSIS */

SATURATION_INCHANNEL = Channel.fromPath("${VEP_SATURATION}/*.vep.gz")

/* CLUSTL */

CLUSTL_OUT = Channel.fromPath("${INTOGEN_DATASETS}/steps/oncodriveclustl/*.clusters_results.tsv")
CLUSTL_OUT.into{ CLUSTL_OUT1; CLUSTL_OUT2 }

process GroupFeaturesCLUSTL {
	tag "Group features OncodriveCLUSTL"
	label "boostdm"
	publishDir "${OUTPUT}/features_group", mode: 'copy'

	input:
        path inputs from CLUSTL_OUT1.collect()
        path cohorts from COHORTS_SUMMARY1

    output:
        path output into CLUSTL_GROUP

	script:
		output = "clustl.tsv.gz"
		"""
		runner.sh features/group.py group-clustl \
			--output ${output} \
			--threshold 0.05 \
			--cohorts ${cohorts} \
			${inputs}
		"""
}

CLUSTL_GROUP_V = CLUSTL_GROUP.first()

OUT_ONCODRIVECLUSTL_CLUSTERS = CLUSTL_OUT2.map{it -> [it.baseName.split('\\.')[0], it]}


/* HOTMAPS */

HOTMAPS_OUT = Channel.fromPath("${INTOGEN_DATASETS}/steps/hotmaps/*.clusters.gz")
HOTMAPS_OUT.into{ HOTMAPS_OUT1; HOTMAPS_OUT2 }

process GroupFeaturesHotMAPS {
	tag "Group features HotMAPS"
	label "boostdm"
	publishDir "${OUTPUT}/features_group", mode: 'copy'

	input:
        path inputs from HOTMAPS_OUT1.collect()
        path cohorts from COHORTS_SUMMARY2

    output:
        path output into HOTMAPS_GROUP

	script:
		output = "hotmaps.tsv.gz"
		"""
		runner.sh features/group.py group-hotmaps \
			--output ${output} \
			--threshold 0.05 \
			--cohorts ${cohorts} \
			${inputs}
		"""
}

HOTMAPS_GROUP_V = HOTMAPS_GROUP.first()

OUT_HOTMAPS_CLUSTERS = HOTMAPS_OUT2.map{it -> [it.baseName.split('\\.')[0], it]}

process RenameHotMAPS {
	tag "Rename HotMAPS clusters file ${cohort}"

	input:
        tuple val(cohort), path(input) from OUT_HOTMAPS_CLUSTERS

    output:
        tuple val(cohort), path(output) into OUT_HOTMAPS_CLUSTERS_RENAMED

	script:
		output = "${cohort}.hotmapsclusters.gz"
		"""
		ln -s ${input} ${output}
		"""
}

/* SMREGIONS */

SMREGIONS_OUT = Channel.fromPath("${INTOGEN_DATASETS}/steps/smregions/*.smregions.tsv.gz")
SMREGIONS_OUT.into{ SMREGIONS_OUT1; SMREGIONS_OUT2 }

process GroupFeaturesSMRegions {
	tag "Group features SMRegions"
	label "boostdm"
	publishDir "${OUTPUT}/features_group", mode: 'copy'

	input:
        path inputs from SMREGIONS_OUT1.collect()
        path cohorts from COHORTS_SUMMARY3

    output:
        path output into SMREGIONS_GROUP

	script:
		output = "smregions.tsv.gz"
		"""
		runner.sh features/group.py group-smregions \
			--output ${output} \
			--threshold 0.1 \
			--cohorts ${cohorts} \
			${inputs}
		"""
}

SMREGIONS_GROUP_V = SMREGIONS_GROUP.first() // convert to value channel


process RenameSMREGIONS {
	tag "Rename SMRegions file ${cohort}"

	input:
        path input from SMREGIONS_OUT2

    output:
        tuple val(cohort), path(output) into OUT_SMREGIONS

	script:
		cohort = input.baseName.split('\\.')[0]
		output = "${cohort}.smregions.gz"
		"""
		ln -s ${input} ${output}
		"""
}

/* DNDS */

DNDS_FILES = Channel.fromPath("${INTOGEN_DATASETS}/steps/dndscv/*.dndscv.tsv.gz")
OUT_DNDSCV = DNDS_FILES.map{it -> [it.baseName.split('\\.')[0], it]}
DNDS_ANNOTMUTS_FILES = Channel.fromPath("${INTOGEN_DATASETS}/steps/dndscv/*.dndscv_annotmuts.tsv.gz")
DNDS_ANNOTMUTS_FILES.into{ DNDS_ANNOTMUTS_FILES1; DNDS_ANNOTMUTS_FILES2 }
OUT_DNDSCV_ANNOTMUTS = DNDS_ANNOTMUTS_FILES1.map{it -> [it.baseName.split('\\.')[0], it]}


/* MutRate */

MUTRATE_FILES = Channel.fromPath("${INTOGEN_DATASETS}/steps/boostDM/mutrate/*.mutrate.json")
OUT_MUTRATE = MUTRATE_FILES.map{it -> [it.baseName.split('\\.')[0], it]}


/* Create datasets */

CREATE_DATASETS_INCHANNEL = OUT_DNDSCV.join(OUT_DNDSCV_ANNOTMUTS).join(OUT_MUTRATE)

process CreateDatasets {
    tag "Creating datasets ${cohort}"
    label "boostdm"
    publishDir "${OUTPUT}/create_datasets", mode: 'copy'

    input:
        tuple val(cohort), path(dndscv), path(dndscvAnnotMuts), path(mutrate) from CREATE_DATASETS_INCHANNEL
        path driversSummary from DRIVERS_SUMMARY_V
        path groupCLUSTL from CLUSTL_GROUP_V
        path groupHotMAPS from HOTMAPS_GROUP_V
        path groupSMRegions from SMREGIONS_GROUP_V

    output:
        tuple val(cohort), path(output) into IN_CV

	script:
		output = "${cohort}.regression_data.tsv"
		"""
		runner.sh annotations/cohort.py \
			--cohort ${cohort} \
			--dndscv-path ${dndscv} \
			--dndscv-annotmuts-path ${dndscvAnnotMuts} \
			--mutrate-path ${mutrate} \
			--clustl-group-path ${groupCLUSTL} \
			--hotmaps-group-path ${groupHotMAPS} \
			--smregions-group-path ${groupSMRegions} \
			--splits ${params.boostdm.bootstrapSplits} \
			--threshold ${params.boostdm.xsThresh} \
                        --out ${output}
		"""
}


process SplitCV {
    tag "Cross validation splits ${cohort}"
    label "boostdm"
    publishDir "${OUTPUT}/splitcv", mode: 'copy'

    input:
        tuple val(cohort), path(input) from IN_CV

    output:
        tuple val(cohort), path(output) into OUT_CV

	script:
		output = "${cohort}.cvdata.pickle.gz"
		// TODO add random seed
		"""
		runner.sh cvdata/cohort.py \
			--input_path ${input} \
			--output_path ${output} \
			--splits ${params.boostdm.bootstrapSplits} \
			--cv ${params.boostdm.cvFraction}
		"""
}

// get only file names
OUT_CV_COHORTS = OUT_CV.map{ it -> it[1] }


process SplitCVMetacohort {
    tag 'Creating the cross validation splits for metacohorts'
    label "boostdm"
    publishDir "${OUTPUT}", mode: 'copy'

    input:
        path(input) from OUT_CV_COHORTS.collect()

    output:
        path("splitcv_meta/*/*.cvdata.pickle.gz") into OUT_CV_META

	script:
		"""
		runner.sh cvdata/meta.py \
			--cores ${task.cpus} \
			--output_path splitcv_meta \
			--input_path .
		"""
}

OUT_CV_META_R = OUT_CV_META.flatten().map{it -> [it.parent.baseName, it.baseName.split('\\.')[0], it]}


/* Training models */

process TrainingMeta {
    tag "Training model ${ttype}-${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/training_meta", mode: 'copy'

    memory { task.memory * task.attempt }
    errorStrategy { task.attempt > 3 ? 'ignore': 'retry' }



	input:
		tuple val(ttype), val(gene), path(input) from OUT_CV_META_R

	output:
		tuple val(ttype), val(gene), path(output) optional true into OUT_TRAIN_META

	script:
		// TODO add seed
		output = "${ttype}/${gene}.models.pickle.gz"
		"""
		mkdir ${ttype}
		runner.sh training.py \
			--splits ${input} \
			--cores ${task.cpus} \
			--min-rows ${params.boostdm.minimumRows} \
			--output ${output}
		"""
}


/* Cross-validation */


process CrossValidation {
    tag "Evaluating model ${ttype}-${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/evaluation", mode: 'copy'

    memory { task.memory * task.attempt }
    errorStrategy { task.attempt > 3 ? 'ignore': 'retry' }


    input:
        tuple val(ttype), val(gene), path(input) from OUT_TRAIN_META

    output:
        tuple val(ttype), val(gene), path(output) into OUT_EVAL_META

	script:
		output = "${ttype}/${gene}.eval.pickle.gz"
		"""
		mkdir ${ttype}
		runner.sh evaluation/auto.py \
			--model ${input} \
			--output ${output}
		"""
}
