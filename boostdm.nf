#!/usr/bin/env nextflow

OUTPUT = "/workspace/projects/boostdm_ch/output/run_20230802_includeLowQuality"

INTOGEN_RUN = "/workspace/datasets/intogen/runs/20200703_oriolRun/CH_IMPACT_out/intogen_merge_20220325"

COHORTS_SUMMARY = Channel.fromPath("${INTOGEN_RUN}/cohorts.tsv")
COHORTS_SUMMARY.into{ COHORTS_SUMMARY1; COHORTS_SUMMARY2; COHORTS_SUMMARY3; COHORTS_SUMMARY4 }

CLUSTL_OUT = Channel.fromPath("${INTOGEN_RUN}/debug/oncodriveclustl/*.tsv")
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
			--ttypes ${cohorts} \
			${inputs}
		"""
}


CLUSTL_GROUP_V = CLUSTL_GROUP.first()


HOTMAPS_OUT = Channel.fromPath("${INTOGEN_RUN}/debug/hotmaps/*.clusters.gz")
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
			--ttypes ${cohorts} \
			${inputs}
		"""
}


HOTMAPS_GROUP_V = HOTMAPS_GROUP.first()


SMREGIONS_OUT = Channel.fromPath("${INTOGEN_RUN}/debug/smregions/*.smregions.tsv.gz")
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
			--ttypes ${cohorts} \
			${inputs}
		"""
}


SMREGIONS_GROUP_V = SMREGIONS_GROUP.first()


DRIVERS_SUMMARY = Channel.fromPath("${INTOGEN_RUN}/drivers.tsv")
DRIVERS_SUMMARY_V = DRIVERS_SUMMARY.first()


DNDS_BASE = "${INTOGEN_RUN}/debug/dndscv"
DNDS_FILES = Channel.fromPath("${DNDS_BASE}/*.dndscv.tsv.gz")
OUT_DNDSCV = DNDS_FILES.map{it -> [it.baseName.split('\\.')[0], it]}


DNDS_ANNOTMUTS_FILES = Channel.fromPath("${DNDS_BASE}/*.dndscv_annotmuts.tsv.gz")
DNDS_ANNOTMUTS_FILES.into{ DNDS_ANNOTMUTS_FILES1; DNDS_ANNOTMUTS_FILES2 }
OUT_DNDSCV_ANNOTMUTS = DNDS_ANNOTMUTS_FILES1.map{it -> [it.baseName.split('\\.')[0], it]}


MUTRATE_FILES = Channel.fromPath("${INTOGEN_RUN}/mutrate/*.json")
OUT_MUTRATE = MUTRATE_FILES.map{it -> [it.baseName.split('\\.')[0], it]}


OUT_ONCODRIVECLUSTL_CLUSTERS = CLUSTL_OUT2.map{it -> [it.baseName.split('\\.')[0], it]}
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


CREATE_DATASETS_INCHANNEL = OUT_DNDSCV.join(OUT_DNDSCV_ANNOTMUTS).join(OUT_MUTRATE)


process CreateDatasets {
    tag "Creating datasets ${cohort}"
    label "boostdm"
    publishDir "${OUTPUT}/create_datasets", mode: 'copy'

    input:
	tuple val(cohort), path(dndscv), path(dndscvAnnotMuts), path(mutrate) from CREATE_DATASETS_INCHANNEL
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
		"""
		runner.sh cvdata/cohort.py \
			--input_path ${input} \
			--output_path ${output} \
			--splits ${params.boostdm.bootstrapSplits} \
			--cv ${params.boostdm.cvFraction}
		"""
}


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
		// TODO set input-path to the list of all files
		"""
		runner.sh cvdata/meta.py \
			--cores ${task.cpus} \
			--output_path splitcv_meta \
			--input_path .
		"""
}


OUT_CV_META_R = OUT_CV_META.flatten().map{it -> [it.parent.baseName, it.baseName.split('\\.')[0], it]}



process TrainingMeta {
    tag "Training model ${ttype}-${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/training_meta", mode: 'copy'

	input:
		tuple val(ttype), val(gene), path(input) from OUT_CV_META_R

	output:
		tuple val(ttype), val(gene), path(output) into OUT_TRAIN_META

	script:
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


process AutoEvaluationMeta {
    tag "Evaluating model ${ttype}-${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/evaluation", mode: 'copy'

    // memory { task.memory * task.attempt }
    // errorStrategy { task.attempt > 3 ? 'ignore': 'retry' }

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


OUT_EVAL_META.into{OUT_EVAL_META1; OUT_EVAL_META2}


process Mutations4Discovery {
    tag 'Mutations for discovery'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        path input from DNDS_ANNOTMUTS_FILES2.collect()

    output:
        path(output) into DISCOVERY_MUTS

	script:
	  	output = "mutations.tsv"
		"""
		runner.sh discovery_index/muts.py \
			--output ${output} \
			${input}
		"""
}


VARIANTS_JSONS_FILES = Channel.from([["boostdm-ch", "${INTOGEN_RUN}/variants.json"]])


process RenameVariants {
	tag "Rename HotMAPS clusters file ${cohort}"

	input:
        tuple val(name), path(input) from VARIANTS_JSONS_FILES

    output:
        path(output) into VARIANTS_JSONS

	script:
		output = "${name}.json"
		"""
		ln -s ${input} ${output}
		"""
}


process Samples4Discovery {
    tag 'Samples for discovery'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        path input from VARIANTS_JSONS.collect()
        path cohorts from COHORTS_SUMMARY4

    output:
        path(output) into DISCOVERY_SAMPLES

	script:
	  	output = "samples.json"
		"""
		runner.sh discovery_index/samples.py \
			--output ${output} \
			--cohorts ${cohorts} \
			${input}
		"""
}


process DiscoveryIndex {
    tag 'Discovery index'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        val (input) from OUT_EVAL_META1.collect()
        path samples from DISCOVERY_SAMPLES
        path mutations from DISCOVERY_MUTS

    output:
        path(output) into DISCOVERY_INDEX

	script:
	  	output = "discovery.tsv"
		"""
		runner.sh discovery_index/discovery.py \
			--output ${output} \
			--mutations ${mutations} \
			--samples ${samples} \
			--evaluation-path ${OUTPUT}/evaluation
		"""
}


process GetEvaluation {
    tag 'Model Selection'
    label "boostdm"
    publishDir "${OUTPUT}/model_selection", mode: 'copy'

    input:
		path discovery from DISCOVERY_INDEX
		val eval from OUT_EVAL_META2.collect()

    output:
        path output into MODEL_OUT

	script:
		output = "eval_data.json"
		"""
		runner.sh evaluation/data.py \
			--eval_folder ${OUTPUT}/evaluation \
			--discovery_path ${discovery} \
			--output ${output}
		"""
}


MODEL_OUT_V = MODEL_OUT.first()

SATURATION = Channel.fromPath("${INTOGEN_RUN}/saturation/maf/*.vep.gz")


process AnnotateSaturation {
    tag "Annotate saturation ${input}"
    label "boostdm"
    publishDir "${OUTPUT}/saturation/annotation", mode: 'copy'

    input:
        path input from SATURATION
        path driversSummary from DRIVERS_SUMMARY_V
        path groupCLUSTL from CLUSTL_GROUP_V
        path groupHotMAPS from HOTMAPS_GROUP_V
        path groupSMRegions from SMREGIONS_GROUP_V

    output:
        tuple val(gene), val(tumor), path(output) into ANNOTATION_SATURATION

    script:
        gene = input.baseName.split('\\.')[0]
        tumor = input.baseName.split('\\.')[1]
        output = "${tumor}/${gene}.annotated.out.gz"
		"""
		mkdir ${tumor}
		runner.sh annotations/gene.py \
			--mutations ${input} \
			--drivers-summary ${driversSummary} \
			--clustl-group ${groupCLUSTL} \
			--hotmaps-group ${groupHotMAPS} \
			--smregions-group ${groupSMRegions} \
			--ttype ${tumor} \
			--output ${output}
		"""
}



process PredictSaturation {
    tag "Predict saturation ${gene} ${ttype}"
    label "boostdm"
    publishDir "${OUTPUT}/saturation/prediction", mode: 'copy'

    input:
        tuple val(gene), val(ttype), path(input) from ANNOTATION_SATURATION
        path model from MODEL_OUT_V

    output:
        path output into PREDICTION_SATURATION

    script:
        output = "${gene}.${ttype}.prediction.tsv.gz"
        """
        runner.sh perform_predictions.py \
                --muts ${input} \
                --gene ${gene} \
                --tumor-type ${ttype} \
                --models-folder ${OUTPUT}/training_meta \
                --evaluations-folder ${OUTPUT}/evaluation \
                --model-selection ${model} \
                --output-file ${output}
        """
}
