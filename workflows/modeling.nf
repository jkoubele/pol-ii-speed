def model_run_id = "model_run_" + new Date().format("yyyy_dd_MM_HH_mm_ss") // Used as name of model results directory

process CreateDesignMatrix {
    input:
    path samplesheet
    val design_formula
    path factor_references

    output:
    path("design_matrix.csv"), emit: design_matrix

    publishDir "${params.outdir}/design_matrix", mode: 'copy'

    script:
    """
    generate_design_matrix.R \
    --samplesheet $samplesheet \
    --formula $design_formula \
     ${factor_references ? "--factor_references $factor_references" : ""}
    """
}

process SplitGeneNames {
    input:
    path gene_names

    output:
    path("gene_names_chunk_*.csv"), emit: gene_names_chunks

    publishDir "${params.outdir}/gene_names_chunks", mode: 'copy'

    script:
    """
    split_gene_names.R \
    --input_gene_names $gene_names \
    --output_folder . \
    --chunk_size 20
    """
}


process RunModel {
    input:
    tuple (
        path(gene_names),
        val(chunk_name),
        path(design_matrix),
        path(exon_counts_matrix),
        path(intron_counts_matrix),
        path(library_size_factors),
        path(coverage_parquet_files)
    )


    output:
    path("model_results_*.csv"), emit: model_result_chunks

    publishDir "${params.outdir}/chunk_model_results/${chunk_name}", mode: 'copy'

    script:
    """
    export PYTHONPATH='${baseDir}'\${PYTHONPATH:+:\$PYTHONPATH}
    run_model.py \
    --design_matrix $design_matrix \
    --gene_names $gene_names \
    --exon_counts $exon_counts_matrix \
    --intron_counts $intron_counts_matrix \
    --library_size_factors $library_size_factors \
    --coverage_data_folder . \
    --output_folder . \
    --output_basename model_results_${chunk_name}
    """

}

process MergeModelResultChunks {
    input:
    path model_result_chunks
    val design_formula

    output:
    path("model_results.csv"), emit: model_results
    path ("design_formula.txt")

    publishDir "${params.outdir}/model_results/${model_run_id}", mode: 'copy'

    script:
    """
    merge_model_result_chunks.R \
    --input_folder . \
    --output_folder . \
    --output_file_name model_results.csv

    echo ${design_formula} > design_formula.txt
    """

}

workflow modeling_workflow {
    take:
        samplesheet
        gene_names_file
        exon_counts
        intron_counts
        library_size_factors
        coverage_files
        design_formula
        factor_reference_levels


    main:
        def samplesheet_channel = Channel.value(file(samplesheet))
        def factor_reference_channel = factor_reference_levels ? Channel.value(file(factor_reference_levels)) : Channel.value([])


        def design_matrix_channel = CreateDesignMatrix(samplesheet_channel,
                                                       design_formula,
                                                       factor_reference_channel)

        gene_names_split = SplitGeneNames(gene_names_file)

        def model_input = gene_names_split.gene_names_chunks
            .flatten()
            .map{ file -> tuple(file, file.baseName)}
            .combine(design_matrix_channel.design_matrix)
            .combine(exon_counts)
            .combine(intron_counts)
            .combine(library_size_factors)
            // Wrapping coverage_files in an extra list prevents unwanted flattening behavior in .combine()
            .combine(coverage_files.map { file_list -> tuple([file_list]) })

        model_result_chunks = RunModel(model_input).model_result_chunks
        collected_results_chunks = model_result_chunks.collect()
        MergeModelResultChunks(collected_results_chunks, design_formula)

}