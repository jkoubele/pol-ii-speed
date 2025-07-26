
process FastQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample), path(read1), path(read2), val(strand) // val(strand) is optional, since we dont use it in this process

    tag "$sample"  // For better logging

    output:
    path("*.zip"), emit: fastqc_reports // we name the output of .zip files so it can be used downstream, e.g. in MultiQC
    path("*.html")

    publishDir "results/fastqc", mode: 'symlink'


    script:
    """
    fastqc $read1 $read2 --outdir .
    """ // Nextflow searches for declared output only in the 'current' folder; since FastQC creates subfolder, we need to move it
}

process MultiQC {
    container 'multiqc/multiqc:latest'

    input:
    path zip_files

    output:
    path("multiqc_report.html")
    path("multiqc_report_data")

    publishDir "results/multiqc", mode: 'copy'

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}



workflow {

//  Channel
//         .fromPath("/cellfile/datapublic/jkoubele/dev-pol-ii-analysis/nextflow_pipeline/results/fastqc/*.zip")
//         .collect()
//         .set { fastqc_zips }

//     MultiQC(fastqc_zips)



    def samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample      = row.sample
            def fq1         = file("${params.fastq_dir}/${row.fastq_1}")
            def fq2         = file("${params.fastq_dir}/${row.fastq_2}")
            def strand      = row.strandedness

            if (!fq1.exists()) error "FASTQ file not found: ${fq1}"
            if (!fq2.exists()) error "FASTQ file not found: ${fq2}"
            if (!['forward', 'reverse'].contains(strand)) {
                error "Invalid strandedness '${strand}' for sample '${sample}'"
            }

            tuple(sample, fq1, fq2, strand)
        }

        def fastqc_out = FastQC(samples)
        fastqc_out.fastqc_reports.collect().set { fastqc_zips }
        MultiQC(fastqc_zips)


}


