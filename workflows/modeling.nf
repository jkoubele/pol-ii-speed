
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

}