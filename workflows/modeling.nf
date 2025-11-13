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
    --formula "$design_formula" \
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
    --chunk_size 100
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
        path(isoform_length_factors),
        path(coverage_parquet_files),
        val(intron_specific_lfc)
    )

    output:
    path("model_results_*.csv"), optional: true, emit: model_result_chunks
    path("logs/ignored_genes.csv"),   emit: ignored_genes_logs
    path("logs/ignored_introns.csv"), emit: ignored_introns_logs

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
    --isoform_length_factors $isoform_length_factors \
    --coverage_data_folder . \
    --intron_specific_lfc ${intron_specific_lfc} \
    --output_folder . \
    --output_basename model_results_${chunk_name}
    """

}

process MergeModelResultChunks {
    input:
    path model_result_chunks
    path ignored_genes_logs
    path ignored_introns_logs
    val design_formula
    val intron_specific_lfc

    output:
    path("model_results.csv"), emit: model_results
    path("model_results_raw.csv")
    path ("model_parameters.txt")
    path("ignored_genes.csv")
    path("ignored_introns.csv")

    publishDir "${params.outdir}/model_results/${model_run_id}", mode: 'copy'

    script:
    """
    merge_and_clean_result_chunks.R \
    --input_folder . \
    --output_folder .

# Intentionally missing indent, so EOF works properly
cat > model_parameters.txt <<EOF
design_formula: "${design_formula}"
intron_specific_lfc: ${intron_specific_lfc}
EOF
    """

}

process CreateVolcanoPlots {
    input:
    path model_results

    output:
    path("volcano_plots/*")

    publishDir "${params.outdir}/model_results/${model_run_id}", mode: 'copy'

    script:
    """
    create_volcano_plots.R \
    --model_results $model_results \
    --output_folder ./volcano_plots \
    """
}

workflow modeling_workflow {
    take:
        samplesheet_input
        gene_names_file_input
        exon_counts_input
        intron_counts_input
        library_size_factors_input
        isoform_length_factors_input
        coverage_files_input
        design_formula
        factor_reference_levels
        intron_specific_lfc


    main:
        def samplesheet_channel = Channel.value(file(samplesheet_input))
        def factor_reference_channel = factor_reference_levels ? Channel.value(file(factor_reference_levels)) : Channel.value([])


        def design_matrix_channel = CreateDesignMatrix(samplesheet_channel,
                                                       design_formula,
                                                       factor_reference_channel)

        gene_names_split = SplitGeneNames(gene_names_file_input)

        def model_input = gene_names_split.gene_names_chunks
            .flatten()
            .map{ file -> tuple(file, file.baseName)}
            .combine(design_matrix_channel.design_matrix)
            .combine(exon_counts_input)
            .combine(intron_counts_input)
            .combine(library_size_factors_input)
            .combine(isoform_length_factors_input)
            // Wrapping coverage_files_input in an extra list prevents unwanted flattening behavior in .combine()
            .combine(coverage_files_input.map { file_list -> tuple([file_list]) })
            // Raw intron_specific_lfc is bool and cannot be used in combine()
            .combine(Channel.value(intron_specific_lfc))

        def run_model_out = RunModel(model_input)

        def collected_results_chunks       = run_model_out.model_result_chunks.filter { it != null }.collect()
        def collected_ignored_genes_logs   = run_model_out.ignored_genes_logs.collect()
        def collected_ignored_introns_logs = run_model_out.ignored_introns_logs.collect()

        def model_result_merged = MergeModelResultChunks(
            collected_results_chunks,
            collected_ignored_genes_logs,
            collected_ignored_introns_logs,
            design_formula,
            intron_specific_lfc
        )

        CreateVolcanoPlots(model_result_merged.model_results)

}