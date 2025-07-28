
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

process ExtractIntronsFromGTF {
    container 'bioinfo_tools'

    input:
    path gtf

    output:
    path("introns.bed"), emit: introns_bed_file

    publishDir "${params.outdir}/extracted_introns", mode: 'copy'


    script:
    """
    extract_introns_from_gtf.py \\
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

process BuildSalmonIndex {
    container 'bioinfo_tools'

    input:
    path transcriptome_fasta

    output:
    path("salmon_index"), emit: salmon_index_dir

    publishDir "${params.outdir}/salmon_index", mode: 'symlink'
    // TODO: Consider using decoy when building the index
    script:
    def gencode_flag = params.gtf_source == 'gencode' ? '--gencode' : ''
    """
    mkdir salmon_index
    salmon index \
      -t $transcriptome_fasta \
      -i salmon_index \
      -p 8 \
      ${gencode_flag}
    """
}

process STARAlign {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(read1), path(read2), val(strand), path(star_index)

    output:
        tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), val(strand), emit: star_bam
        path("${sample}.Log.final.out")

    publishDir "${params.outdir}/star/${sample}", mode: 'copy'

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

process ExtractIntronicReads {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(bam_file), val(strand), path(introns_bed_file)

    output:
        tuple val(sample), path("intronic_reads_sorted.bam"), emit:  intronic_bams
        tuple val(sample), path("intronic_reads_plus_strand.bed.gz"), path("intronic_reads_minus_strand.bed.gz"),  emit:  intronic_bed_files
        tuple val(sample), path("intron_read_counts.tsv"), emit:  intron_read_counts

    publishDir "${params.outdir}/intronic_reads/${sample}", mode: 'copy'

    tag "$sample"

    script:
    """
    extract_intronic_reads.py \\
        --input_bam $bam_file \\
        --intron_bed_file $introns_bed_file \\
        --strandedness $strand

    sort -k 1,1 -k 2,2n intronic_reads_plus_strand.bed > tmp_plus_strand.bed
    mv tmp_plus_strand.bed intronic_reads_plus_strand.bed

    sort -k 1,1 -k 2,2n intronic_reads_minus_strand.bed > tmp_minus_strand.bed
    mv tmp_minus_strand.bed intronic_reads_minus_strand.bed

    pigz -f intronic_reads_plus_strand.bed
    pigz -f intronic_reads_minus_strand.bed
    """
}

process RemoveIntronicReadsFromFASTQ {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(read1), path(read2), path(bam_introns)

    output:
        tuple val(sample), path("R1.fastq.gz"), path("R2.fastq.gz"), emit: exonic_fastq

    publishDir "${params.outdir}/FASTQ_without_intronic_reads/${sample}", mode: 'copy'

    tag "$sample"

    script:
    """
    remove_intronic_reads_from_fastq.py \\
        --bam_introns $bam_introns \\
        --input_fastq $read1 \\
        --output_fastq R1.fastq
    pigz -f R1.fastq

    remove_intronic_reads_from_fastq.py \\
        --bam_introns $bam_introns \\
        --input_fastq $read2 \\
        --output_fastq R2.fastq
    pigz -f R2.fastq
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

    def salmon_index_channel
    if (params.salmon_index) {
        log.info "Using provided Salmon index: ${params.salmon_index}"
        salmon_index_channel = Channel.value(file(params.salmon_index))
    } else {
        log.info "No Salmon index provided. Building from transcriptome FASTA."
        salmon_index_channel = BuildSalmonIndex(file(params.transcriptome_fasta)).salmon_index_dir
    }

    def introns = ExtractIntronsFromGTF(file(params.gtf_file))

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
    fastqc_out.fastqc_reports.collect() | MultiQC

    aligned_bams = samples
    .combine(star_index_channel)
    .map { sample, r1, r2, strand, index -> tuple(sample, r1, r2, strand, index) }
    | STARAlign


    extracted_intronic_reads = aligned_bams.star_bam
    .combine(introns.introns_bed_file)
    .map { sample, bam, strand, bed -> tuple(sample, bam, strand, bed) }
    | ExtractIntronicReads

    exonic_fastq = samples
    .map { sample, r1, r2, strand -> tuple(sample, [r1, r2]) }
    .join(extracted_intronic_reads.intronic_bams).map { sample, reads, intronic_bam ->
        tuple(sample, reads[0], reads[1], intronic_bam)
    } | RemoveIntronicReadsFromFASTQ

//     exonic_fastq.exonic_fastq | view

}
