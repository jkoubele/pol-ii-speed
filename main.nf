include { preprocessing_workflow } from './workflows/preprocessing.nf'
include { modeling_workflow } from './workflows/modeling.nf'

workflow{
    if (!['preprocess', 'model', 'all'].contains(params.stage)){
     error "Invalid 'stage' value: ${params.stage}. Must be either 'preprocess', 'model' or 'all'."
    }

    if (['preprocess', 'all'].contains(params.stage)) {
        // Run pre-processing workflow
        if (!params.samplesheet ||
        !params.fastq_dir ||
        !params.gtf_file ||
        !params.genome_fasta ||
        !params.transcriptome_fasta ||
        params.salmon_index_with_decoy == null ||
        !params.gtf_source) {
            error "Missing required parameters! Please run with -params-file <your_params.yaml>."
        }

        if (!['ensembl', 'gencode'].contains(params.gtf_source)) {
            error "Invalid 'gtf_source' value: ${params.gtf_source}. Must be 'ensembl' or 'gencode'."
        }

        preprocessing_workflow(
            params.samplesheet,
            params.fastq_dir,
            params.gtf_file,
            params.genome_fasta,
            params.transcriptome_fasta,
            params.gtf_source,
            params.salmon_index_with_decoy,
            params.star_index,
            params.salmon_index)

    }

    if (['model', 'all'].contains(params.stage)) {
        modeling_workflow(
        params.samplesheet,
        null,
        null,
        null,
        null,
        null,
        params.design_formula,
        params.factor_reference_levels)
    }
}





