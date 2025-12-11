library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(pbapply)
library(BiocParallel)


pboptions(type = "txt") 
register(MulticoreParam(workers = 4, progressbar = FALSE))
gtf <- rtracklayer::import("/cellfile/datapublic/jkoubele/reference_genomes/ensembl_115_GRCm39/Mus_musculus.GRCm39.115.gtf")

exons <- gtf[gtf$type == "exon"]
# exons <- exons[exons$transcript_biotype == "protein_coding"]

# Split by gene
exons_by_gene <- split(exons, exons$gene_id)

t0 <- Sys.time()
constitutive_exons_by_gene <- bplapply(exons_by_gene[1:100], function(gene_exons) {
  transcripts_ids <- unique(gene_exons$transcript_id)
  num_transcripts <- length(transcripts_ids)
  
  atomic_exon_ranges <- disjoin(gene_exons, with.revmap = TRUE, ignore.strand = FALSE)
  
  atomic_exon_to_transcripts <- mcols(atomic_exon_ranges)$revmap
  
  keep_mask <- vapply(atomic_exon_to_transcripts, function(exon_transcript_indices) {
    length(unique(gene_exons$transcript_id[exon_transcript_indices])) == num_transcripts
  }, logical(1))
  atomic_exon_ranges[keep_mask]
})

t1 <- Sys.time()
print(t1 - t0)

