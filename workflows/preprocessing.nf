def preprocessing_output_subfolder = "preprocessing"

process FastQC {
    input:
    tuple val(sample), path(read1), path(read2)

    tag "$sample"

    output:
    path("*.zip"), emit: fastqc_reports
    path("*.html")

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/FastQC", mode: 'copy'

    script:
    """
    fastqc $read1 $read2 --outdir .
    """
}

process MultiQC {
    input:
    path zip_files

    output:
    path("multiqc_report.html")
    path("multiqc_report_data")

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/MultiQC", mode: 'copy'

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}


process ExtractGenomicFeatures {
    input:
    path gtf

    output:
        path('introns.bed'), emit: introns_bed_file
        path('constitutive_exons.bed'), emit: constitutive_exons_bed_file
        path('introns.gtf'), emit: introns_gtf_file
        path('constitutive_exons.gtf'), emit: constitutive_exons_gtf_file
        path("protein_coding_genes.csv"), emit: protein_coding_gene_names
        path("all_genes.csv")

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/genomic_features", mode: 'copy'

    script:
    """
    extract_genomic_features_from_gtf.R \
        --gtf $gtf \
        --threads ${task.cpus}
    """
}


process BuildStarIndex {
    input:
    path fasta
    path gtf

    output:
    path("star_index"), emit: star_index_dir

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/STAR_index", mode: 'copy'

    script:
    """
    mkdir star_index
    STAR \
      --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir star_index \
      --genomeFastaFiles $fasta \
      --sjdbGTFfile $gtf \
      --limitGenomeGenerateRAM 190000000000
    """
}

process BuildSalmonIndex {
    input:
    path transcriptome_fasta
    path genome_fasta
    path gtf
    val(gtf_source)
    val(salmon_index_with_decoy)

    output:
    path("salmon_index"), emit: salmon_index_dir
    path("decoy_transcriptome/gentrome.fa", optional: true)
    path("decoy_transcriptome/decoys.txt", optional: true)

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/Salmon_index", mode: 'copy'

    script:
    def gencode_flag = gtf_source == 'gencode' ? '--gencode' : ''
    if (salmon_index_with_decoy) {
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
    } else {
        """
        mkdir salmon_index
        salmon index \
            -t $transcriptome_fasta \
            -i salmon_index \
            -p ${task.cpus} \
            ${gencode_flag}
        """
    }
}

process PrepareTx2Gene {
    input:
    path gtf

    output:
    path("tx2gene.tsv"), emit: tx2gene_file

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/tx2gene", mode: 'copy'

    script:
    """
    prepare_tx2gene.R --gtf_file $gtf
    """
}

process CreateGenomeFastaIndex {
    input:
    path genome_fasta

    output:
    path("*.fai"), emit: genome_fai_file

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/genome_fai", mode: 'copy'

    script:
    """
    samtools faidx $genome_fasta
    """
}

process STARAlign {
    input:
        tuple val(sample), path(read1), path(read2), path(star_index)

    output:
        tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: star_bam
        path("${sample}.*")

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/STAR/${sample}", mode: 'copy'

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
      --limitBAMsortRAM 60000000000
    """
}

process ExtractIntronicReads {
    input:
        tuple val(sample), path(bam_file), val(strand), path(introns_bed_file)

    output:
        tuple val(sample), path("intronic_reads_sorted.bam"), emit:  intronic_bam_files
        tuple val(sample), path("intronic_reads_plus_strand.bed.gz"), path("intronic_reads_minus_strand.bed.gz"),  emit:  intronic_bed_files
        tuple val(sample), path("${sample}.intron_read_counts.tsv"), emit:  intron_read_counts

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/intronic_reads/${sample}", mode: 'copy'

    tag "$sample"

    script:
    """
    extract_intronic_reads.py \\
        --input_bam $bam_file \\
        --intron_bed_file $introns_bed_file \\
        --strandedness $strand

    mv intron_read_counts.tsv ${sample}.intron_read_counts.tsv

    sort -k 1,1 -k 2,2n intronic_reads_plus_strand.bed > tmp_plus_strand.bed
    mv tmp_plus_strand.bed intronic_reads_plus_strand.bed

    sort -k 1,1 -k 2,2n intronic_reads_minus_strand.bed > tmp_minus_strand.bed
    mv tmp_minus_strand.bed intronic_reads_minus_strand.bed

    pigz -p ${task.cpus} -f intronic_reads_plus_strand.bed
    pigz -p ${task.cpus} -f intronic_reads_minus_strand.bed
    """
}

process SalmonQuantification {
    input:
        tuple val(sample), path(read1), path(read2), path(salmon_index)

    output:
        tuple val(sample), path("${sample}.quant.sf"), emit: salmon_quant
        path("**")

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/Salmon_quantification/${sample}", mode: 'copy'

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
    mv quant.sf ${sample}.quant.sf
    """
}

process ConstitutiveExonsFeatureCounts {
    input:
        tuple val(sample), path(bam), val(strand), path(constitutive_exons_gtf_file)

    output:
        tuple val(sample), path("${sample}.constitutive_exons.txt"), emit: constitutive_exon_counts
        path "${sample}.constitutive_exons.txt.summary"

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/constitutive_exons_read_count/${sample}", mode: 'copy'

    tag "$sample"

    script:
    def strandedness_argument = (strand == 'forward') ? 1 : (strand == 'reverse') ? 2 : 0
    """
    featureCounts \
      -a ${constitutive_exons_gtf_file} \
      -t constitutive_exon \
      -g constitutive_exon_id \
      -p -s ${strandedness_argument} -B -C -T ${task.cpus} \
      -o ${sample}.constitutive_exons.txt \
      ${bam}

    # Rename last header column (named after the BAM file by default) to 'NumReads'
    awk 'BEGIN{OFS="\t"} \$0 !~ /^#/ && !done { \$NF="NumReads"; done=1 } { print }' \
     ${sample}.constitutive_exons.txt > ${sample}.constitutive_exons.tmp \
     && mv ${sample}.constitutive_exons.tmp ${sample}.constitutive_exons.txt
    """
}

process ComputeCoverage {
    input:
        tuple val(sample), path(bed_file_plus), path(bed_file_minus), path(genome_fai_file)

    output:
        tuple val(sample), path("coverage_plus.bedGraph.gz"), path("coverage_minus.bedGraph.gz"), emit: bed_graph_files

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/bed_graph_intron_coverage/${sample}", mode: 'copy'

    tag "$sample"

    script:
    """
    bedtools genomecov \
        -bga \
        -split \
        -i $bed_file_plus \
        -g $genome_fai_file \
        > coverage_plus.bedGraph
        pigz -p ${task.cpus} -f coverage_plus.bedGraph

    bedtools genomecov \
        -bga \
        -split \
        -i $bed_file_minus \
        -g $genome_fai_file \
        > coverage_minus.bedGraph
        pigz -p ${task.cpus} -f coverage_minus.bedGraph
    """
}

process RescaleCoverage {
    input:
        tuple val(sample), path(bedgraph_file_plus), path(bedgraph_file_minus), path(introns_bed_file)


    output:
    path("${sample}.parquet"), emit: coverage_parquet_file

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/rescaled_coverage", mode: 'copy'

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

process AggregateReadCounts {
    input:
    tuple val(sample_names), path(exon_quant_files), path(intron_counts_files), path(tx2gene)

    output:
        path("exon_counts.tsv"), emit: exon_counts
        path("intron_counts.tsv"), emit: intron_counts
        path("library_size_factors.tsv"), emit: library_size_factors
        path("isoform_length_factors.tsv"), emit: isoform_length_factors

    publishDir "${params.outdir}/${preprocessing_output_subfolder}/aggregated_counts", mode: 'copy'

    script:
    """
    aggregate_read_counts.R \
      --tx2gene $tx2gene \
      --sample_names ${sample_names.join(' ')} \
      --exon_quant_files ${exon_quant_files.join(' ')} \
      --intron_counts_files ${intron_counts_files.join(' ')}
    """
}


workflow preprocessing_workflow {
    take:
        samplesheet
        fastq_dir
        gtf_file
        genome_fasta
        transcriptome_fasta
        gtf_source
        salmon_index_with_decoy
        star_index
        salmon_index

    main:
        def gtf_channel = Channel.value(file(gtf_file))
        def genome_fasta_channel = Channel.value(file(genome_fasta))
        def transcriptome_fasta_channel = Channel.value(file(transcriptome_fasta))

        def star_index_channel
        if (star_index) {
            log.info "Using provided STAR index: ${star_index}"
            star_index_channel = Channel.value(file(star_index))
        } else {
            log.info "No STAR index provided. Building from genome FASTA and GTF."
            star_index_channel = BuildStarIndex(genome_fasta_channel, gtf_channel).star_index_dir
        }

        def salmon_index_channel
        if (salmon_index) {
            log.info "Using provided Salmon index: ${salmon_index}"
            salmon_index_channel = Channel.value(file(salmon_index))
        } else {
            log.info "No Salmon index provided. Building from transcriptome FASTA."
            salmon_index_channel = BuildSalmonIndex(
                transcriptome_fasta_channel,
                genome_fasta_channel,
                gtf_channel,
                gtf_source,
                salmon_index_with_decoy).salmon_index_dir
        }

        samples = Channel
            .fromPath(samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def sample      = row.sample
                def fq1         = file("${fastq_dir}/${row.fastq_1}")
                def fq2         = file("${fastq_dir}/${row.fastq_2}")
                def strand      = row.strandedness

                if (!fq1.exists()) error "FASTQ file not found: ${fq1}"
                if (!fq2.exists()) error "FASTQ file not found: ${fq2}"
                if (!['forward', 'reverse'].contains(strand)) {
                    error "Invalid strandedness '${strand}' for sample '${sample}'"
                }

                tuple(sample, fq1, fq2, strand)
            }


        def tx2gene_out = PrepareTx2Gene(gtf_channel)
        def fai_index = CreateGenomeFastaIndex(genome_fasta_channel).genome_fai_file

        def genomic_features = ExtractGenomicFeatures(gtf_channel)

        def fastqc_out = samples.map{sample, fq1, fq2, strand -> tuple(sample, fq1, fq2)} | FastQC
        def fastqc_out_aggregated = fastqc_out.fastqc_reports.collect()
        fastqc_out_aggregated | MultiQC

        def aligned_bams = samples
        .combine(star_index_channel)
        .map { sample, r1, r2, strand, star_index ->
            tuple(sample, r1, r2, star_index)
        } | STARAlign

        def constitutive_exons_read_counts = aligned_bams.star_bam
       .join(samples.map { sample, r1, r2, strand -> tuple(sample, strand) })
       .combine(genomic_features.constitutive_exons_gtf_file)
       | ConstitutiveExonsFeatureCounts

       def extracted_intronic_reads = aligned_bams.star_bam
       .join(samples.map { sample, r1, r2, strand -> tuple(sample, strand) })
       .combine(genomic_features.introns_bed_file)
       | ExtractIntronicReads

       def salmon_quant_out = samples
       .map { sample, r1, r2, strand -> tuple(sample, r1, r2) }
       .combine(salmon_index_channel) | SalmonQuantification

       def bed_graph_coverage = extracted_intronic_reads.intronic_bed_files
       .combine(fai_index)| ComputeCoverage

       def rescaled_coverage = bed_graph_coverage.bed_graph_files
       .combine(genomic_features.introns_bed_file)| RescaleCoverage

       def rescaled_coverage_combined = rescaled_coverage.collect()

       def data_aggregation =  salmon_quant_out.salmon_quant
       .join(extracted_intronic_reads.intron_read_counts)
       .collect(flat: false).map { list_of_tuples ->
            def list_of_tuples_sorted = list_of_tuples.sort { it[0] }
            def sample_names = list_of_tuples_sorted*.getAt(0)
            def quant_files  = list_of_tuples_sorted*.getAt(1)
            def intron_files = list_of_tuples_sorted*.getAt(2)
            tuple(sample_names, quant_files, intron_files)
       }
       .combine(tx2gene_out.tx2gene_file) | AggregateReadCounts

    emit:
        gene_names_file          = genomic_features.protein_coding_gene_names
        exon_counts              = data_aggregation.exon_counts
        intron_counts            = data_aggregation.intron_counts
        library_size_factors     = data_aggregation.library_size_factors
        isoform_length_factors   = data_aggregation.isoform_length_factors
        coverage_files           = rescaled_coverage_combined

}