include { preprocessing_workflow } from './workflows/preprocessing.nf'


def path_aggregated_counts = "${params.outdir}/aggregated_counts"
def path_rescaled_coverage = "${params.outdir}/rescaled_coverage"
def path_gene_names = "${params.outdir}/gene_names"


workflow{

    if !['preprocess', 'model', 'all'].contains(params.stage){
     error "Invalid 'stage' value: ${params.stage}. Must be either 'preprocess', 'model' or 'all'."
    }

    if ['preprocess', 'all'].contains(params.stage) {
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
            samplesheet              : params.samplesheet,
            fastq_dir                : params.fastq_dir,
            gtf_file                 : params.gtf_file,
            genome_fasta             : params.genome_fasta,
            transcriptome_fasta      : params.transcriptome_fasta,
            gtf_source               : params.gtf_source,
            salmon_index_with_decoy  : params.salmon_index_with_decoy,
            star_index               : params.star_index,
            salmon_index             : params.salmon_index
            path_aggregated_counts   : path_aggregated_counts
            path_rescaled_coverage   : path_rescaled_coverage
            path_gene_names          : path_gene_names
        )
    }

    if ['model', 'all'].contains(params.stage) {

    }


}
