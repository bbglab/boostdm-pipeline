#!/usr/bin/env nextflow


// Path channels

DNDS_ANNOTMUTS_FILES = Channel.fromPath("${INTOGEN_DATASETS}/steps/dndscv/*.dndscv_annotmuts.tsv.gz")
COHORTS_SUMMARY = Channel.fromPath("${INTOGEN_DATASETS}/cohorts.tsv")
OUT_EVAL_PATH = Channel.fromPath("${OUTPUT}/evaluation/*/*.eval.pickle.gz")
OUT_EVAL = OUT_EVAL_PATH.map{it -> [it.getParent().baseName, it.baseName.split('\\.')[0], it]}


process Mutations4Discovery {

    tag 'Mutations for discovery'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        path input from DNDS_ANNOTMUTS_FILES.collect()

    output:
        path(output) into DISCOVERY_MUTS

	script:
	  	output = "mutations.tsv.gz"
		"""
		runner.sh discovery_index/muts.py \
			--output ${output} \
			${input}
		"""
}


VARIANTS_STATS_JSON_FOLDER = Channel.fromPath("${INTOGEN_DATASETS}/steps/variants", type: 'dir')


process CollectVariants {

    tag 'Create variants.json'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        path input from VARIANTS_STATS_JSON_FOLDER

    output:
        path(output) into VARIANTS_JSON

    script:
  	output = "variants.json"
	"""
	runner.sh discovery_index/preprocess_variants.py \
		--inputfolder ${input} \
		--output ${output}
	"""
}


process Samples4Discovery {

    tag 'Samples for discovery'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        path input from VARIANTS_JSON
        path cohorts from COHORTS_SUMMARY

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
        val (input) from OUT_EVAL.collect()
        path samples from DISCOVERY_SAMPLES
        path mutations from DISCOVERY_MUTS

    output:
        path(output) into DISCOVERY_INDEX

	script:
	  	output = "discovery.tsv.gz"
		"""
		runner.sh discovery_index/discovery.py \
			--output ${output} \
			--mutations ${mutations} \
			--samples ${samples} \
			--evaluation-path ${OUTPUT}/evaluation
		"""
}
