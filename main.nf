include { preprocessing_workflow } from './workflows/preprocessing.nf'
include { modeling_workflow } from './workflows/modeling.nf'


workflow {

    def preproc_out

    if (!['preprocess', 'model', 'all'].contains(params.stage)) {
        error("Invalid 'stage' value: ${params.stage}. Must be either 'preprocess', 'model' or 'all'.")
    }

    if (['preprocess', 'all'].contains(params.stage)) {
        // Run pre-processing workflow
        if (!params.samplesheet || !params.fastq_dir || !params.gtf_file || !params.genome_fasta || !params.transcriptome_fasta || params.salmon_index_with_decoy == null || !params.gtf_source) {
            error("Missing required parameters! Please run with -params-file <your_params.yaml>.")
        }

        if (!['ensembl', 'gencode'].contains(params.gtf_source)) {
            error("Invalid 'gtf_source' value: ${params.gtf_source}. Must be 'ensembl' or 'gencode'.")
        }

        preproc_out = preprocessing_workflow(
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
        // Run modeling workflow
        def gene_names_file
        def exon_counts
        def intron_counts
        def library_size_factors
        def isoform_length_factors
        def coverage_files

        if (params.stage == 'model') {
            def preprocessing_output_subfolder = "preprocessing"
            gene_names_file        = Channel.value(file("${params.outdir}/${preprocessing_output_subfolder}/gene_names/protein_coding_genes.csv"))
            exon_counts            = Channel.value(file("${params.outdir}/${preprocessing_output_subfolder}/aggregated_counts/exon_counts.tsv"))
            intron_counts          = Channel.value(file("${params.outdir}/${preprocessing_output_subfolder}/aggregated_counts/intron_counts.tsv"))
            library_size_factors   = Channel.value(file("${params.outdir}/${preprocessing_output_subfolder}/aggregated_counts/library_size_factors.tsv"))
            isoform_length_factors = Channel.value(file("${params.outdir}/${preprocessing_output_subfolder}/aggregated_counts/isoform_length_factors.tsv"))
            coverage_files         = Channel.fromPath("${params.outdir}/${preprocessing_output_subfolder}/rescaled_coverage/*.parquet").collect()
        }

        if (params.stage == 'all') {

            gene_names_file        = preproc_out.gene_names_file
            exon_counts            = preproc_out.exon_counts
            intron_counts          = preproc_out.intron_counts
            library_size_factors   = preproc_out.library_size_factors
            isoform_length_factors = preproc_out.isoform_length_factors
            coverage_files         = preproc_out.coverage_files
        }


        def lrt_contrasts = params.lrt_contrasts ?: []

        modeling_workflow(
            params.samplesheet,
            gene_names_file,
            exon_counts,
            intron_counts,
            library_size_factors,
            isoform_length_factors,
            coverage_files,
            params.design_formula,
            lrt_contrasts,
            params.intron_specific_lfc
        )
    }
}
