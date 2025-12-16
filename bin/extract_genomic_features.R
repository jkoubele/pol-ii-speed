#!/usr/bin/env Rscript

library(argparse)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)
library(readr)

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


args <- parser$parse_args()

register(SnowParam(workers = args$threads, type = "SOCK", progressbar = FALSE))

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

gtf <- rtracklayer::import(args$gtf)

genes <- gtf[gtf$type == "gene"]
exons <- gtf[gtf$type == "exon"]

genes <- genes[!is.na(genes$gene_id)]
exons <- exons[!is.na(exons$gene_id)]

# Extract gene IDs
gene_type_column <- if ("gene_biotype" %in% names(mcols(genes))) {
  "gene_biotype"
} else if ("gene_type" %in% names(mcols(genes))) {
  "gene_type"
} else {
  stop("Neither 'gene_biotype' nor 'gene_type' found in GTF attributes.")
}

all_gene_ids <- sort(unique(genes$gene_id))
protein_coding_gene_ids <- sort(unique(genes$gene_id[mcols(genes)[[gene_type_column]] == "protein_coding"]))

readr::write_csv(data.frame(gene_id = all_gene_ids),
                 file.path(output_folder, "all_genes.csv"))
readr::write_csv(data.frame(gene_id = protein_coding_gene_ids),
                 file.path(output_folder, "protein_coding_genes.csv"))

# Different annotations (Ensemble vs. Gencode) use different naming for UTRs
utr_types <- intersect(c("five_prime_utr", "three_prime_utr", "UTR"), unique(gtf$type))

exons_and_utr <- gtf[gtf$type %in% c("exon", utr_types)]
exons_and_utr <- exons_and_utr[!is.na(exons_and_utr$gene_id)]

exons_by_gene <- split(exons, exons$gene_id)
gene_ranges_by_gene <- split(genes, genes$gene_id)
exons_and_utr_by_gene <- split(exons_and_utr, exons_and_utr$gene_id)

# Extract introns
introns_by_gene <- bplapply(
  names(gene_ranges_by_gene),
  function(gene_id, gene_ranges_by_gene, exons_and_utr_by_gene) {

    gene_range <- range(gene_ranges_by_gene[[gene_id]])
    exons_and_utr_ranges <- exons_and_utr_by_gene[[gene_id]]

    if (is.null(exons_and_utr_ranges) || length(exons_and_utr_ranges) == 0) {
      warning("No exon/UTR for gene ", gene_id)
      return(GRanges())
    }

    exons_and_utr_ranges <- GenomicRanges::reduce(exons_and_utr_ranges, ignore.strand = FALSE)
    GenomicRanges::setdiff(gene_range, exons_and_utr_ranges, ignore.strand = FALSE)
  },
  gene_ranges_by_gene = gene_ranges_by_gene,
  exons_and_utr_by_gene = exons_and_utr_by_gene
)
names(introns_by_gene) <- names(gene_ranges_by_gene)


# Extract constitutive exons
constitutive_exons_by_gene <- bplapply(exons_by_gene, function(gene_exons) {
  if (length(gene_exons) == 0) return(GRanges())

  transcripts_ids <- unique(gene_exons$transcript_id)
  num_transcripts <- length(transcripts_ids)

  atomic_exon_ranges <- disjoin(gene_exons, with.revmap = TRUE, ignore.strand = FALSE)

  atomic_exon_to_transcripts <- mcols(atomic_exon_ranges)$revmap

  keep_mask <- vapply(atomic_exon_to_transcripts, function(exon_transcript_indices) {
    length(unique(gene_exons$transcript_id[exon_transcript_indices])) == num_transcripts
  }, logical(1))
  atomic_exon_ranges[keep_mask]
})


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

  # Assign number in the 5' to 3' direction
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
