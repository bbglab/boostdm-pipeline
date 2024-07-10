SATURATION_OUT = Channel.fromPath("${OUTPUT}/saturation/prediction/*.model.*.features.*.prediction.tsv.gz")
GENE_TTYPE_SATURATION_OUT = SATURATION_OUT.map{ it -> [it.baseName.split('\\.')[0], it.baseName.split('\\.')[2], it.baseName.split('\\.')[4]]}


process Blueprints {
    tag "Generate blueprints for gene=${gene} model=${ttype_model} features=${ttype_features}"
    label "boostdm"
    publishDir "${OUTPUT}/output_plots/blueprints", mode: 'copy'

    input:
        tuple val(gene), val(ttype_model), val(ttype_features) from GENE_TTYPE_SATURATION_OUT
        
    output:
        path("${gene}.model.${ttype_model}.features.${ttype_features}.*") into OUT_BLUEPRINT

    script:
        
        """
        runner.sh output_plots/blueprint.py \
                  --gene ${gene} \
                  --ttmodel ${ttype_model} \
                  --ttfeatures ${ttype_features}
        """
}


process ClusteredBlueprints {
    tag "Generate all blueprint stacks with hierarchical clustering"
    label "boostdm"
    publishDir "${OUTPUT}/output_plots/clustered_blueprints", mode: 'copy'

    output:
        path("*.clustered_blueprint.*") into OUT_CLUSTERED_BLUEPRINT

    script:
        
        """
        runner.sh output_plots/clustered_blueprint.py
        """
}


TRAINING_OUT = Channel.fromPath("${OUTPUT}/training_meta/*/*.models.pickle.gz")
GENE_TTYPE_TRAINING_OUT = TRAINING_OUT.map{ it -> [it.baseName.split('\\.')[0], it.getParent().baseName]}

process DiscoveryBending {
    tag "Generate inverse exponential fit with unique mutation count in subsamples for gene=${gene} ttype=${ttype}"
    label "boostdm"
    publishDir "${OUTPUT}/output_plots/discovery_bending", mode: 'copy'

    input:
        tuple val(gene), val(ttype) from GENE_TTYPE_TRAINING_OUT

    output:
        path("*.*.bending.*") into OUT_DISCOVERY_BENDING

    script:
        
        """
        runner.sh output_plots/discovery_plot.py \
                  --gene ${gene} \
                  --ttmodel ${ttype}
        """
}
