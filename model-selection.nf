/* Model selection */


EVAL_OUT = Channel.fromPath("${OUTPUT}/evaluation/*/*.eval.pickle.gz")
OUT_EVAL_META = EVAL_OUT.collect()


process ModelSelection {
    tag 'Model Selection -- which models have to be run for which gene-ttypes'
    label "boostdm"
    publishDir "${OUTPUT}/model_selection", mode: 'copy'

    input:
        val(input) from OUT_EVAL_META

    output:
        path(output) into MODEL

        script:
                output = "eval_data.pickle.gz"
                """
                runner.sh evaluation/data.py \
                        --eval_folder ${OUTPUT}/evaluation \
                        --discovery_path ${OUTPUT}/discovery/discovery.tsv.gz \
                        --output ${output}
                """
}

