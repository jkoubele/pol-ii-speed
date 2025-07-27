
process FastQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample), path(read1), path(read2), val(strand) // val(strand) is optional, since we dont use it in this process

    tag "$sample"

    output:
    path("*.zip"), emit: fastqc_reports
    path("*.html")

    publishDir "${params.outdir}/fastqc", mode: 'symlink'


    script:
    """
    fastqc $read1 $read2 --outdir .
    """
}

process MultiQC {
    container 'multiqc/multiqc:latest'

    input:
    path zip_files

    output:
    path("multiqc_report.html")
    path("multiqc_report_data")

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

process ExtractIntrons {
    container 'bioinfo_tools'

    input:
    path gtf

    output:
    path("introns.bed"), emit: introns_bed_file

    publishDir "${params.outdir}/extracted_introns", mode: 'copy'


    script:
    """
    python3 ${projectDir}/scripts/extract_introns_from_gtf.py \\
        --gtf_file $gtf \\
        --gtf_source ${params.gtf_source} \\
    """
}

process BuildStarIndex {
    container 'bioinfo_tools'

    input:
    path fasta
    path gtf

    output:
    path("star_index"), emit: star_index_dir

    publishDir "${params.outdir}/STAR_index", mode: 'symlink'

    script:
    """
    mkdir star_index
    STAR \
      --runThreadN 16 \
      --runMode genomeGenerate \
      --genomeDir star_index \
      --genomeFastaFiles $fasta \
      --sjdbGTFfile $gtf
    """
}

process STARAlign {
    container 'bioinfo_tools'

    input:
    tuple val(sample), path(read1), path(read2), val(strand), path(star_index)

    output:
    path("${sample}.Aligned.sortedByCoord.out.bam"), emit: star_bam
    path("${sample}.Log.final.out")

    publishDir "${params.outdir}/star", mode: 'copy'

    tag "$sample"

    script:
    """
    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir $star_index \
      --readFilesIn $read1 $read2 \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes All \
      --quantMode GeneCounts \
      --peOverlapNbasesMin 1 \
      --outFileNamePrefix ${sample}. \
      --limitBAMsortRAM 50000000000
    """
}



workflow {
    def star_index_channel
    if (params.star_index) {
        log.info "Using provided STAR index: ${params.star_index}"
        star_index_channel = Channel.value(file(params.star_index))
    } else {
        log.info "No STAR index provided. Building from genome FASTA and GTF."
        star_index_channel = BuildStarIndex(file(params.genome_fasta), file(params.gtf_file)).star_index_dir
    }

    def introns = ExtractIntrons(file(params.gtf_file))

    samples = Channel
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

    fastqc_out = FastQC(samples)
    fastqc_out.fastqc_reports.collect().set { fastqc_zips }
    MultiQC(fastqc_zips)

    samples.combine(star_index_channel) | STARAlign

}
