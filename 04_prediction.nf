#!/usr/bin/env nextflow


CLUSTL_GROUP = Channel.fromPath("${OUTPUT}/features_group/clustl.tsv.gz")
CLUSTL_GROUP_V = CLUSTL_GROUP.first()

HOTMAPS_GROUP = Channel.fromPath("${OUTPUT}/features_group/hotmaps.tsv.gz")
HOTMAPS_GROUP_V = HOTMAPS_GROUP.first()

SMREGIONS_GROUP = Channel.fromPath("${OUTPUT}/features_group/smregions.tsv.gz")
SMREGIONS_GROUP_V = SMREGIONS_GROUP.first()

SPLITCV_OUT = Channel.fromPath("${OUTPUT}/splitcv_meta/*/*.models.pickle.gz")
GENE_TTYPE_OUT = SPLITCV_OUT.map{ it -> [it.baseName.split('\\.')[0], it.getParent().baseName]}


process AnnotateSaturation {
    tag "Annotate saturation ${gene} ${ttype}"
    label "boostdm"
    publishDir "${OUTPUT}/saturation/annotation", mode: 'copy'

    input:
        tuple val(gene), val(ttype) from GENE_TTYPE_OUT
        path groupCLUSTL from CLUSTL_GROUP_V
        path groupHotMAPS from HOTMAPS_GROUP_V
        path groupSMRegions from SMREGIONS_GROUP_V

    output:
        tuple val(gene), val(ttype), path(output) into ANNOTATION_SATURATION

    script:
        vep = "${VEP_SATURATION}/${gene}.vep.gz"
        output = "${gene}.${ttype}.annotated.tsv.gz"
        
        """
        runner.sh annotations/gene.py \
                --gene ${gene} \
                --ttype ${ttype} \
                --mutations ${vep} \
                --clustl-group ${groupCLUSTL} \
                --hotmaps-group ${groupHotMAPS} \
                --smregions-group ${groupSMRegions} \
        """
}


MODEL = Channel.fromPath("${OUTPUT}/model_selection/eval_data.pickle.gz").first()

// this is a toy change

// for testing purposes only
// SATURATION_OUT = Channel.fromPath("${OUTPUT}/saturation/annotation/*.annotated.tsv.gz")
// ANNOTATION_SATURATION = SATURATION_OUT.map{ it -> [it.baseName.split('\\.')[0], it.baseName.split('\\.')[1], it]}


process PredictSaturation {
    tag "Predict saturation ${gene} ${ttype}"
    label "boostdm"
    publishDir "${OUTPUT}/saturation/prediction", mode: 'copy'

    input:
        tuple val(gene), val(ttype), path(input) from ANNOTATION_SATURATION
        path model from MODEL

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
                --output-file ${output} \
                --high-quality-only
        """
}


