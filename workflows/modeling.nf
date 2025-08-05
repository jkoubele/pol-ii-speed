
process CreateDesignMatrix {
    container 'pol_ii_bioconductor'

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
    container 'pol_ii_bioconductor'

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
    --chunk_size 1000
    """
}


process RunModel {
    container 'bioinfo_tools'

    input:
    path design_matrix
    path gene_names
    path exon_counts_matrix
    path intron_counts_matrix
    path library_size_factors
    path coverage_parquet_files

    output:
    path("success.txt"), emit: success_message

    publishDir "${params.outdir}/model_results", mode: 'copy'

    script:
    """
    test_load_dataset.py \
    --design_matrix $design_matrix \
    --gene_names $gene_names \
    --exon_counts $exon_counts_matrix \
    --intron_counts $intron_counts_matrix \
    --library_size_factors $library_size_factors \
    --coverage_data_folder .
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
        def design_matrix_channel = CreateDesignMatrix(samplesheet_channel, design_formula, factor_reference_channel)

        gene_names_split = SplitGeneNames(gene_names_file)
        gene_names_split.gene_names_chunks | view


        def model_results = RunModel(
            design_matrix_channel.design_matrix,
            gene_names_file,
            exon_counts,
            intron_counts,
            library_size_factors,
            coverage_files
            )
}