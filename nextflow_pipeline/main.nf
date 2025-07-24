
process FastQC {
    tag "$sample"  // For better logging

    input:
    tuple val(sample), path(read1), path(read2), val(strand)

    output:
    path("*.zip"), emit: fastqc_reports
    path("*.html")

    script:
    """
    mkdir -p fastqc_output
    fastqc $read1 $read2 --outdir fastqc_output
    mv fastqc_output/* .
    """
}


workflow {
    def samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample      = row.sample
            def fq1         = file("${params.fastq_dir}/${row.fastq_1}")
            def fq2         = file("${params.fastq_dir}/${row.fastq_2}")
            def strand      = row.strandedness

            // Validate input
            if (!fq1.exists()) error "FASTQ file not found: ${fq1}"
            if (!fq2.exists()) error "FASTQ file not found: ${fq2}"

            if (!['forward', 'reverse'].contains(strand)) {
                error "Invalid strandedness '${strand}' for sample '${sample}'. Only 'forward' or 'reverse' allowed."
            }

            tuple(sample, fq1, fq2, strand)
        }


    // TEST process: print sample info
    samples | view //{ it -> "Sample: ${it[0]}, R1: ${it[1].name}, R2: ${it[2].name}, Strand: ${it[3]}" }
}
