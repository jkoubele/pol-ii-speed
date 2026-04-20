def modeling_output_subfolder = "modeling"
def model_run_id = "model_run_" + new Date().format("dd_MMM_yyyy_HH_mm_ss")


process WriteDesignMetadata {
    input:
    val design_formula
    val intron_specific_lfc
    val lrt_contrasts

    output:
    path "design.json"
    path "lrt_contrasts.json", emit: lrt_contrasts_json

    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/design_metadata", mode: 'copy'

    script:
    def design_obj = [
        design_formula      : design_formula,
        intron_specific_lfc : intron_specific_lfc
    ]

    def design_json_str       = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(design_obj))
    def lrt_contrasts_json_str = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(lrt_contrasts ?: []))

    """
    printf '%s\n' '${design_json_str}'        > design.json
    printf '%s\n' '${lrt_contrasts_json_str}' > lrt_contrasts.json
    """
}



process CreateDesignMatrices {
    input:
    path samplesheet
    val design_formula
    path lrt_contrasts_json

    output:
    path("design_matrix.tsv"), emit: design_matrix
    path("lrt_metadata.tsv"), emit: lrt_metadata
    path("reduced_design_matrices"), emit: reduced_design_matrices

    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/design_matrices", mode: 'copy'

    script:
    """
    create_design_matrices.R \
    --samplesheet $samplesheet \
    --formula "$design_formula" \
    --lrt_contrasts_json $lrt_contrasts_json
    """
}

process GetModeledGenes {
    input:
    tuple(path(gene_names), path(modelable_genes), path(modelable_introns))

    output:
    path("gene_chunks/modeled_genes_chunk_*.tsv"), emit: gene_names_chunks
    path("modeled_genes.tsv")
    path("modeled_introns.tsv"), emit: modeled_introns

    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/modeled_genes", mode: 'copy'

    script:
    """
    get_modeled_genes.R \
    --input_gene_names $gene_names \
    --modelable_genes $modelable_genes \
    --modelable_introns $modelable_introns \
    --output_folder . \
    --chunk_size 100
    """
}


process FitModel {
    input:
    tuple (
        path(modeled_genes_chunk),
        val(chunk_name),
        path(modeled_introns),
        path(design_matrix),
        path(lrt_metadata),
        path(reduced_matrices_folder),
        path(exon_counts_matrix),
        path(intron_counts_matrix),
        path(library_size_factors),
        path(isoform_length_factors),
        path(coverage_parquet_files),
        val(intron_specific_lfc)
    )

    output:
    path("model_parameters*.tsv"), emit: model_parameters_chunk
    path("test_results*.tsv"), emit: test_results_chunk
    tuple path("cache_for_regularization*.pt"), val(chunk_name), emit: cache_for_regularization

//     publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/processing_chunks/model_results_chunks/${chunk_name}", mode: 'copy'

    script:
    """
    export PYTHONPATH='${baseDir}'\${PYTHONPATH:+:\$PYTHONPATH}
    fit_model.py \
    --modeled_genes $modeled_genes_chunk \
    --modeled_introns $modeled_introns \
    --design_matrix $design_matrix \
    --exon_counts $exon_counts_matrix \
    --intron_counts $intron_counts_matrix \
    --library_size_factors $library_size_factors \
    --isoform_length_factors $isoform_length_factors \
    --coverage_data_folder . \
    --lrt_metadata ${lrt_metadata} \
    --reduced_matrices_folder ${reduced_matrices_folder} \
    --intron_specific_lfc ${intron_specific_lfc} \
    --output_folder . \
    --output_name_suffix _${chunk_name}
    """

}

process MergeModelResultChunks {
    input:
    path model_parameters_chunks
    path test_results_chunks

    output:
    path("test_results_before_regularization.tsv"), emit: test_results_before_regularization
    path("model_parameters.tsv"), emit: model_parameters
    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/model_results",
     mode: 'copy',
     saveAs: { filename ->
            if( filename == 'test_results_before_regularization.tsv' ) null // Not publishing intermediate test results
            else filename
        }

    script:
    """
    merge_result_chunks.R \
    --input_folder . \
    --output_folder .
    """

}


process AdaptiveShrinkage {
    input:
    path model_parameters

    output:
    path("regularization_coefficients.tsv"), emit: regularization_coefficients

    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/adaptive_shrinkage", mode: 'copy'

    script:
    """
    adaptive_shrinkage.R \
    --model_parameters $model_parameters
    """
}

process FitRegularizedModel {
    input:
    tuple(
        path(cache_for_regularization),
        val(chunk_name),
        path(regularization_coefficients)
    )

    output:
    path("regularized_model_parameters*.tsv"), emit: regularized_model_parameters_chunk

//     publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/processing_chunks/regularized_model_results_chunks/${chunk_name}", mode: 'copy'

    script:
    """
    export PYTHONPATH='${baseDir}'\${PYTHONPATH:+:\$PYTHONPATH}
    fit_regularized_model.py \
    --cache_for_regularization $cache_for_regularization \
    --regularization_coefficients $regularization_coefficients \
    --output_name_suffix _$chunk_name
    """
}

process AddRegularizationToTestResults {
    input:
    path test_results_before_regularization
    path regularized_model_parameters_chunks

    output:
    path("test_results.tsv"), emit: test_results
    path("regularized_model_parameters.tsv")

    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/model_results", mode: 'copy'

    script:
    """
    add_regularization_to_test_results.R \
    --regularization_chunks_folder . \
    --test_results $test_results_before_regularization \
    --output_folder .
    """

}

process CreateVolcanoPlots {
    input:
    path test_results

    output:
    path("**")

    publishDir "${params.outdir}/${modeling_output_subfolder}/${model_run_id}/volcano_plots", mode: 'copy'

    script:
    """
    create_volcano_plots.R \
    --test_results $test_results
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
        lrt_contrasts
        intron_specific_lfc
        modelable_genes_input
        modelable_introns_input


    main:
        def samplesheet = Channel.value(file(samplesheet_input))

        def design_metadata = WriteDesignMetadata(
            design_formula,
            intron_specific_lfc,
            lrt_contrasts
        )

        def design_matrices_data = CreateDesignMatrices(
            samplesheet,
            design_formula,
            design_metadata.lrt_contrasts_json
        )
        get_modeled_genes_input = gene_names_file_input
        .combine(modelable_genes_input)
        .combine(modelable_introns_input)

        def modeled_genes = GetModeledGenes(get_modeled_genes_input)

        def model_input = modeled_genes.gene_names_chunks
            .flatten()
            .map{ file -> tuple(file, file.baseName)}
            .combine(modeled_genes.modeled_introns)
            .combine(design_matrices_data.design_matrix)
            .combine(design_matrices_data.lrt_metadata)
            .combine(design_matrices_data.reduced_design_matrices)
            .combine(exon_counts_input)
            .combine(intron_counts_input)
            .combine(library_size_factors_input)
            .combine(isoform_length_factors_input)
            // Wrapping coverage_files_input in an extra list prevents unwanted flattening behavior in .combine()
            .combine(coverage_files_input.map { file_list -> tuple([file_list]) })
            // Raw intron_specific_lfc is bool and cannot be used in combine()
            .combine(Channel.value(intron_specific_lfc))

        def fit_model_output = FitModel(model_input)
//
//         def collected_model_parameters_chunks  = fit_model_output.model_parameters_chunk.filter { it != null }.collect()
//         def collected_test_results_chunks  = fit_model_output.test_results_chunk.filter { it != null }.collect()
//
//         def model_result_merged = MergeModelResultChunks(
//             collected_model_parameters_chunks,
//             collected_test_results_chunks
//         )
//
//         def adaptive_shrinkage_out = AdaptiveShrinkage(model_result_merged.model_parameters)
//
//         def regularized_model_input = fit_model_output.cache_for_regularization
//             .combine(adaptive_shrinkage_out.regularization_coefficients)
//
//         def fit_regularized_model_output = FitRegularizedModel(regularized_model_input)
//
//         def collected_regularized_model_parameters_chunk = fit_regularized_model_output.regularized_model_parameters_chunk.collect()
//
//         def add_regularization_output = AddRegularizationToTestResults(
//             model_result_merged.test_results_before_regularization,
//             collected_regularized_model_parameters_chunk
//         )
//
//         CreateVolcanoPlots(add_regularization_output.test_results)

}
