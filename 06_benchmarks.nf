// prepare vep input -> benchmarks/vep_input/

process PrepareVEPInput {

    tag "Generate mutation files in input format compatible with VEP"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/vep_input/", mode: 'copy'
        
    output:
        path("*.tsv") into VEP_INPUT

    script:
        
        """
        runner.sh benchmarks/prepare_vep_input.py
        
        """
}

// create cv tables -> benchmarks/cv_tables/

VEP_INPUT.into{ VEP_INPUT_1; VEP_INPUT_2; VEP_INPUT_3 }

VEP_INPUT_TUPLES_1 = VEP_INPUT_1.flatten().map{ it -> [it.baseName.split('\\.')[0], it]}

process CreateCVTables {

    tag "Generate cross-validation testing data for gene=${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/cv_tables/", mode: 'copy'

    input:
        tuple val(gene), path(vep_input) from VEP_INPUT_TUPLES_1
        
    output:
        path("*.*.50.iter.tsv") optional true into OUT_CV_TABLES

    script:
        
        """
        runner.sh benchmarks/create_cv_tables.py --input ${vep_input}

        """
}


// run vep with dbNSPF plugin -> benchmarks/vep_output_dbNSFP

VEP_INPUT_TUPLES_2 = VEP_INPUT_2.flatten().map{ it -> [it.baseName.split('\\.')[0], it]}

process RunVEPdbNSFP{

    tag "Run VEP on saturation with dbNSFP plugin for gene=${gene}"
    label "vep"
    publishDir "${OUTPUT}/benchmarks/vep_output_dbNSFP/", mode: 'copy'
    memory = 20.GB

    input:
        tuple val(gene), path(vep_input) from VEP_INPUT_TUPLES_2
        
    output:
        path(output) into VEP_DBNSFP_OUTPUT

    script:
        output = "${gene}.vep.dbnsfp.tsv"
        vep_dir = "/workspace/datasets/vep/"
        dbnsfp = "/workspace/datasets/vep/homo_sapiens/plugins/dbNSFP4.5a_grch38.gz"
        scores = "SIFT_score,SIFT4G_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationAssessor_score,FATHMM_score,MetaLR_score,MetaRNN_score,CADD_raw,VEST4_score,PROVEAN_score,REVEL_score,ESM1b_score,EVE_score,AlphaMissense_score,phyloP100way_vertebrate,phyloP470way_mammalian,phyloP17way_primate"
    
        """
        vep --dir ${vep_dir} --format ensembl --tab -i ${vep_input} \
            --offline --cache -o ${output} \
            --species homo_sapiens --assembly GRCh38 --fork 8 --mane_select --plugin dbNSFP,${dbnsfp},${scores}

        """

}


VEP_DBNSFP_OUTPUT_TUPLE = VEP_DBNSFP_OUTPUT.map{ it -> [it.baseName.split('\\.')[0], it]}

process ReformatVEP{

    tag "Reformat output from VEP"
    label "boostdm"

    input:
        tuple val(gene), path(input) from VEP_DBNSFP_OUTPUT_TUPLE
        
    output:
        tuple val(gene), path(output) into VEP_DBNSFP_OUTPUT_REFORMATTED

    script:

        output = "${gene}.tsv"

        """
        cat ${input} | grep -v "^##" > ${output}

        """

}


process Saturation_dbNSFP {

    tag "Generate saturation data with dbNSFP annotations for gene=${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/saturation_dbNSFP/", mode: 'copy'

    input:
        tuple val(gene), path(vep_input) from VEP_DBNSFP_OUTPUT_REFORMATTED
        
    output:
        path("*.*.saturation.dbNSFP.tsv") optional true into SATURATION_DBNSFP

    script:
    
        """
        runner.sh benchmarks/saturation_dbnsfp.py --input ${vep_input}

        """

}


VEP_INPUT_TUPLES_3 = VEP_INPUT_3.flatten().map{ it -> [it.baseName.split('\\.')[0], it]}

process Saturation_MAVE {

    tag "Generate saturation data with MAVE annotations for gene=${gene}"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/saturation_mave/", mode: 'copy'

    input:
        tuple val(gene), path(vep_input) from VEP_INPUT_TUPLES_3
        
    output:
        path("*.*.saturation.mave.tsv") optional true into SATURATION_MAVE

    script:
    
        """
        runner.sh benchmarks/saturation_mave.py --input ${vep_input}

        """

}


OUT_CV_TABLES_READY = OUT_CV_TABLES.collect()
SATURATION_DBNSFP_READY = SATURATION_DBNSFP.collect()
SATURATION_MAVE_READY = SATURATION_MAVE.collect()
CV_DBNSFP_READY = OUT_CV_TABLES_READY.concat(SATURATION_DBNSFP_READY, SATURATION_MAVE_READY)


process AnnotateCVTables {

    tag "Annotate CV tables with dbNSFP and MAVE data"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/cv_tables_annotated/", mode: 'copy'

    input:
        val(input) from CV_DBNSFP_READY.collect()
        
    output:
        path("*.*.50.iter.annotated.tsv") into CV_ANNOTATED_OUT

    script:
    
        """
        runner.sh benchmarks/annotate_cv_tables.py

        """

}


CV_ANNOTATED_OUT_READY = CV_ANNOTATED_OUT.collect()
CV_ANNOTATED_OUT_READY.into{ CV_ANNOTATED_OUT_READY_1; CV_ANNOTATED_OUT_READY_2 }


process BenchmarkTable {
    
    tag "Compute all PR-AUC for benchmark data"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/", mode: 'copy'

    input:
        val(input) from CV_ANNOTATED_OUT_READY_1
        
    output:
        path("benchmark.tsv") into BENCHMARK_TABLE

    script:
    
        """
        runner.sh benchmarks/precision_recall.py

        """
}


process BenchmarkPlots {
    
    tag "Dump all PR plots"
    label "boostdm"
    publishDir "${OUTPUT}/benchmarks/pr_plots", mode: 'copy'

    input:
        val(input) from CV_ANNOTATED_OUT_READY_2
        
    output:
        path("*.*.*.prauc.*") into PR_PLOTS

    script:
    
        """
        runner.sh benchmarks/precision_recall.py --plots

        """
}
