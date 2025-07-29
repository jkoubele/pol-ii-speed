
process FastQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample), path(read1), path(read2)

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
        --gtf_source ${params.gtf_source}
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
      --runThreadN ${task.cpus} \
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
    path genome_fasta
    path gtf

    output:
    path("salmon_index"), emit: salmon_index_dir
    path ("decoy_transcriptome/gentrome.fa")
    path("decoy_transcriptome/decoys.txt")

    publishDir "${params.outdir}/salmon_index", mode: 'symlink'

    script:
    def gencode_flag = params.gtf_source == 'gencode' ? '--gencode' : ''
    """
    generateDecoyTranscriptome.sh \
        -a $gtf \
        -g $genome_fasta \
        -t $transcriptome_fasta \
        -o decoy_transcriptome \
        -j ${task.cpus}

    mkdir salmon_index
    salmon index \
      -t decoy_transcriptome/gentrome.fa \
      -d decoy_transcriptome/decoys.txt \
      -i salmon_index \
      -p ${task.cpus} \
      ${gencode_flag}
    """
}

process PrepareTx2Gene {
    container 'pol_ii_bioconductor'

    input:
    path gtf

    output:
    path("tx2gene.tsv"), emit: tx2gene_file

    publishDir "${params.outdir}/tx2gene", mode: 'copy'

    script:
    """
    prepare_tx2gene.R --gtf_file $gtf
    """
}

process CreateGenomeFastaIndex {
    container 'bioinfo_tools'

    input:
    path genome_fasta

    output:
    path("*.fai"), emit: genome_fai_file

    publishDir "${params.outdir}/genome_fai", mode: 'copy'

    script:
    """
    samtools faidx $genome_fasta
    """
}

process STARAlign {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(read1), path(read2), path(star_index)

    output:
        tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: star_bam
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
      --peOverlapNbasesMin 10 \
      --outFileNamePrefix ${sample}. \
      --limitBAMsortRAM 30000000000
    """
}

process ExtractIntronicReads {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(bam_file), val(strand), path(introns_bed_file)

    output:
        tuple val(sample), path("intronic_reads_sorted.bam"), emit:  intronic_bam_files
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

process SalmonQuantification {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(read1), path(read2), path(salmon_index)

    output:
        tuple val(sample), path("quant.sf"), emit: salmon_quant
        tuple path("lib_format_counts.json"), path("logs/salmon_quant.log")

    publishDir "${params.outdir}/salmon_quantification/${sample}", mode: 'copy'

    tag "$sample"

    script:
    """
    salmon quant \
        -i $salmon_index \
        -l A \
        -1 $read1 \
        -2 $read2 \
        -p ${task.cpus} \
        -o .
    """
}

process ComputeCoverage {
    container 'bioinfo_tools'

    input:
        tuple val(sample), path(bed_file_plus), path(bed_file_minus), path(genome_fai_file)

    output:
        tuple val(sample), path("coverage_plus.bedGraph.gz"), path("coverage_minus.bedGraph.gz"), emit: bed_graph_files

    publishDir "${params.outdir}/bed_graph_intron_coverage/${sample}", mode: 'copy'

    tag "$sample"

    script:
    """
    bedtools genomecov \
        -bga \
        -split \
        -i $bed_file_plus \
        -g $genome_fai_file \
        > coverage_plus.bedGraph
        pigz -f coverage_plus.bedGraph

    bedtools genomecov \
        -bga \
        -split \
        -i $bed_file_minus \
        -g $genome_fai_file \
        > coverage_minus.bedGraph
        pigz -f coverage_minus.bedGraph
    """
}

process RescaleCoverage {
    container 'pol_ii_bioconductor'

    input:
        tuple val(sample), path(bedgraph_file_plus), path(bedgraph_file_minus), path(introns_bed_file)

    output:
        tuple val(sample), path("${sample}.parquet"), emit: coverage_parquet_file

    publishDir "${params.outdir}/rescaled_coverage", mode: 'copy'

    tag "$sample"

    script:
    """
    compute_rescaled_coverage.R \
        --bed_graph_plus $bedgraph_file_plus \
        --bed_graph_minus $bedgraph_file_minus \
        --introns_bed_file $introns_bed_file \
        --output_file_basename $sample
    """
}

workflow {

    def gtf_channel = Channel.value(file(params.gtf_file))
    def genome_fasta_channel = Channel.value(file(params.genome_fasta))
    def transcriptome_fasta_channel = Channel.value(file(params.transcriptome_fasta))

    def star_index_channel
    if (params.star_index) {
        log.info "Using provided STAR index: ${params.star_index}"
        star_index_channel = Channel.value(file(params.star_index))
    } else {
        log.info "No STAR index provided. Building from genome FASTA and GTF."
        star_index_channel = BuildStarIndex(genome_fasta_channel, gtf_channel).star_index_dir
    }

    def salmon_index_channel
    if (params.salmon_index) {
        log.info "Using provided Salmon index: ${params.salmon_index}"
        salmon_index_channel = Channel.value(file(params.salmon_index))
    } else {
        log.info "No Salmon index provided. Building from transcriptome FASTA."
        salmon_index_channel = BuildSalmonIndex(transcriptome_fasta_channel, genome_fasta_channel, gtf_channel).salmon_index_dir
    }

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

    def introns_bed_channel = ExtractIntronsFromGTF(gtf_channel).introns_bed_file
    def tx2gene_out = PrepareTx2Gene(gtf_channel).tx2gene_file
    def fai_index = CreateGenomeFastaIndex(genome_fasta_channel).genome_fai_file

    def fastqc_out = samples.map{sample, fq1, fq2, strand -> tuple(sample, fq1, fq2)} | FastQC
    fastqc_out.fastqc_reports.collect() | MultiQC

    def aligned_bams = samples
    .combine(star_index_channel)
    .map { sample, r1, r2, strand, star_index ->
        tuple(sample, r1, r2, star_index)
    } | STARAlign

   def extracted_intronic_reads = aligned_bams.star_bam
   .join(samples.map { sample, r1, r2, strand -> tuple(sample, strand) })
   .combine(introns_bed_channel)
   | ExtractIntronicReads

   def exonic_fastq = samples
   .map { sample, r1, r2, strand -> tuple(sample, r1, r2) }
   .join(extracted_intronic_reads.intronic_bam_files)
   | RemoveIntronicReadsFromFASTQ

   def salmon_quant_out = exonic_fastq.exonic_fastq
   .combine(salmon_index_channel) | SalmonQuantification

   def bed_graph_coverage = extracted_intronic_reads.intronic_bed_files
   .combine(fai_index)| ComputeCoverage

   def rescaled_coverage = bed_graph_coverage.bed_graph_files
   .combine(introns_bed_channel)| RescaleCoverage

}
