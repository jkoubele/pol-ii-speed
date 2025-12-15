library(argparse)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)

parser <- ArgumentParser()
parser$add_argument("--gtf",
                    type = "character",
                    required = TRUE,
                    help = "Path to the input GTF file.")
parser$add_argument("--threads",
                    type = "integer",
                    default = 4,
                    help = "Number of threads used by BiocParallel.")
parser$add_argument("--output_folder",
                    type = "character",
                    default = '.')
parser$add_argument("--min_intron_length",
                    type = "integer",
                    default = 50)
parser$add_argument("--min_constitutive_exon_length",
                    type = "integer",
                    default = 20)

if (interactive()) {
  args <- list(gtf = "/cellfile/datapublic/jkoubele/reference_genomes/ensembl_115_GRCm39/Mus_musculus.GRCm39.115.gtf",
               threads = 4,
               output_folder = "/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/test_genome_features",
               min_intron_length = 50,
               min_constitutive_exon_length = 20)

} else {
  args <- parser$parse_args()
}

register(MulticoreParam(workers = args$threads, progressbar = FALSE))

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

gtf <- rtracklayer::import(args$gtf)

genes <- gtf[gtf$type == "gene"]
exons <- gtf[gtf$type == "exon"]

# Different anotations (Ensemble vs. Gencode) use different naming for UTRs
utr_types <- intersect(c("five_prime_utr", "three_prime_utr", "UTR"), unique(gtf$type))

exons_and_utr <- gtf[gtf$type %in% c("exon", utr_types)]

# Split by gene
exons_by_gene <- split(exons, exons$gene_id)
gene_ranges_by_gene <- split(genes, genes$gene_id)
exons_and_utr_by_gene <- split(exons_and_utr, exons_and_utr$gene_id)

t0 <- Sys.time()
# Extract introns
introns_by_gene <- bplapply(names(gene_ranges_by_gene)[1:1000], function(gene_id) {
  gene_range <- range(gene_ranges_by_gene[[gene_id]])
  exons_and_utr_ranges <- exons_and_utr_by_gene[[gene_id]]

  if (is.null(exons_and_utr_ranges) || length(exons_and_utr_ranges) == 0) {
    warning("No exon/UTR for gene ", gene_id)
    return(GRanges())
  }

  exons_and_utr_ranges <- GenomicRanges::reduce(exons_and_utr_ranges, ignore.strand = FALSE)
  introns <- GenomicRanges::setdiff(gene_range, exons_and_utr_ranges, ignore.strand = FALSE)
  introns
})
names(introns_by_gene) <- names(gene_ranges_by_gene)[1:1000]


# Extract constitutive exons
constitutive_exons_by_gene <- bplapply(exons_by_gene[1:1000], function(gene_exons) {
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

constitutive_exons <- unlist(GRangesList(constitutive_exons_by_gene), use.names = TRUE)
constitutive_exons$gene_id <- names(constitutive_exons)
names(constitutive_exons) <- NULL
constitutive_exons <- constitutive_exons[width(constitutive_exons) >= args$min_constitutive_exon_length]

introns <- unlist(GRangesList(introns_by_gene), use.names = TRUE)
introns$gene_id <- names(introns)
names(introns) <- NULL
introns <- introns[width(introns) >= args$min_intron_length]

features_list <- list(
  introns = introns,
  constitutive_exons = constitutive_exons
)

for (feature_name in names(features_list)) {
  features <- features_list[[feature_name]]
  num_overlapping_genes <- countOverlaps(features, genes, ignore.strand = FALSE)
  features <- features[num_overlapping_genes == 1]

  # Assign exon number in the 5' to 3' direction
  features <- features[order(features$gene_id, start(features), end(features))]
  feature_index_in_gene <- ave(seq_along(features),
                               features$gene_id,
                               FUN = seq_along)
  num_features_in_gene <- ave(seq_along(features),
                              features$gene_id,
                              FUN = length)
  mask_minus_strand <- as.character(strand(features)) == "-"

  features$number_in_gene <- feature_index_in_gene
  features$number_in_gene[mask_minus_strand] <- num_features_in_gene[mask_minus_strand] + 1 - feature_index_in_gene[mask_minus_strand]
  features$name <- paste0(features$gene_id, "_", features$number_in_gene)

  features$score <- 0
  features <- sort(features)
  features_list[[feature_name]] <- features
}

rtracklayer::export(features_list[["introns"]], file.path(output_folder, "introns.bed"), format = "BED")
rtracklayer::export(features_list[["constitutive_exons"]], file.path(output_folder, "constitutive_exons.bed"), format = "BED")
