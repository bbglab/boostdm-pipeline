#!/usr/bin/env nextflow


// Path channels

OUT_EVAL_PATH = Channel.fromPath("${OUTPUT}/evaluation/*/*.eval.pickle.gz")
OUT_EVAL = OUT_EVAL_PATH.map{it -> [it.getParent().baseName, it.baseName.split('\\.')[0], it]}


process DiscoveryIndex {

    tag 'Discovery index'
    label "boostdm"
    publishDir "${OUTPUT}/discovery", mode: 'copy'

    input:
        val (input) from OUT_EVAL.collect()
        
    output:
        path(output) into DISCOVERY_INDEX

	script:
	  	output = "discovery.tsv.gz"
		"""
		runner.sh discovery_index/discovery.py \
			--output ${output} \
			--mutations ${MUTATIONS} \
			--evaluation-path ${OUTPUT}/evaluation
		"""
}
