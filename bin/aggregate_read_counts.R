#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(tximport)
library(DESeq2)


parser <- ArgumentParser()
parser$add_argument("--tx2gene",
                    required = TRUE,
                    help = "Path to tx2gene.tsv")
parser$add_argument("--exon_quant_files",
                    nargs = "+", 
                    required = TRUE)
parser$add_argument("--intron_counts_files", 
                    nargs = "+", 
                    required = TRUE)
parser$add_argument("--constitutive_exon_counts_files", 
                    nargs = "+", 
                    required = TRUE)
parser$add_argument("--sample_names", nargs = "+", required = TRUE)
parser$add_argument("--output_folder", default = '.')

args <- parser$parse_args()

output_folder <- args$output_folder
sample_names <- args$sample_names
exon_quant_files <- args$exon_quant_files
intron_counts_files <- args$intron_counts_files
constitutive_exon_counts_files <- args$constitutive_exon_counts_files

tx2gene <- read_tsv(args$tx2gene, show_col_types = FALSE) |>
  mutate(
    TXNAME = sub("\\|.*$", "", TXNAME),
    TXNAME = sub("\\.[0-9]+$", "", TXNAME)
  )

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

stopifnot(length(sample_names) == length(exon_quant_files),
          length(sample_names) == length(intron_counts_files),
          length(sample_names) == length(constitutive_exon_counts_files))

names(exon_quant_files) <- args$sample_names
txi <- tximport(exon_quant_files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE,
                countsFromAbundance = 'no')

exon_read_counts <- as.data.frame(txi$counts) |>
  rownames_to_column(var = "gene_id")

avg_isoform_length <- as.data.frame(txi$length) |>
  rownames_to_column(var = "gene_id")
write_tsv(avg_isoform_length, file.path(output_folder, "avg_isoform_length.tsv"))

intron_counts_by_sample <- map2(sample_names, intron_counts_files, function(sample, file) {
  read_tsv(file) |>
    select(name, count) |>
    set_names(c("name", sample))
})
intron_read_counts <- purrr::reduce(intron_counts_by_sample, full_join, by = "name") |>
  dplyr::rename(intron_id = name)

constitutive_exon_counts_by_sample <- map2(sample_names, constitutive_exon_counts_files, function(sample, file) {
  read_tsv(file, comment = "#", show_col_types = FALSE) |>
    select(Geneid, NumReads) |>
    set_names(c("constitutive_exon_id", sample))
})

constitutive_exon_read_counts <- purrr::reduce(constitutive_exon_counts_by_sample, full_join, by = "constitutive_exon_id") |>
  dplyr::select(constitutive_exon_id, all_of(sample_names))

exon_read_counts <- exon_read_counts |> select(gene_id, all_of(sample_names))
intron_read_counts <- intron_read_counts |> select(intron_id, all_of(sample_names))

write_tsv(exon_read_counts, file.path(args$output_folder, 'exon_counts.tsv'))
write_tsv(intron_read_counts, file.path(args$output_folder, 'intron_counts.tsv'))
write_tsv(constitutive_exon_read_counts, file.path(args$output_folder, "constitutive_exon_counts.tsv"))

exon_counts_matrix <- exon_read_counts |>
  column_to_rownames("gene_id") |>
  as.matrix()
intron_counts_matrix <- intron_read_counts |>
  column_to_rownames("intron_id") |>
  as.matrix()
all_counts_matrix <- rbind(exon_counts_matrix, intron_counts_matrix) |> round()

avg_isoform_length_matrix <- avg_isoform_length |>
  column_to_rownames("gene_id") |>
  as.matrix()
stopifnot(identical(colnames(avg_isoform_length_matrix), sample_names))

length_introns <- matrix(
  1.0,
  nrow = nrow(intron_counts_matrix),
  ncol = length(sample_names),
  dimnames = list(rownames(intron_counts_matrix), sample_names)
)

avg_length_all <- rbind(
  avg_isoform_length_matrix[rownames(exon_counts_matrix), sample_names, drop = FALSE],
  length_introns
)
stopifnot(identical(rownames(avg_length_all), rownames(all_counts_matrix)))
stopifnot(identical(colnames(avg_length_all), colnames(all_counts_matrix)))

dds <- DESeqDataSetFromMatrix(
  countData = all_counts_matrix,
  colData = DataFrame(row.names = colnames(all_counts_matrix)),
  design = ~1
)
assays(dds)[["avgTxLength"]] <- avg_length_all

dds <- estimateSizeFactors(dds)
normalization_factors <- normalizationFactors(dds)

first_intron_row <- rownames(intron_counts_matrix)[1]
library_size_factors <- normalization_factors[first_intron_row,]
length_size_factors <- exp(log(avg_length_all) - rowMeans(log(avg_length_all)))
stopifnot(identical(colnames(length_size_factors), names(library_size_factors)))

reconstructed_normalization_factors <- t(t(length_size_factors) * library_size_factors)
stopifnot(max(abs(normalization_factors - reconstructed_normalization_factors)) < 1e-8)

length_size_factors_genes <- length_size_factors[rownames(exon_counts_matrix),]

library_size_factors |>
  enframe(name = "sample", value = "library_size_factor") |>
  write_tsv(file.path(output_folder, "library_size_factors.tsv"))

length_size_factors_genes |>
  as.data.frame() |>
  rownames_to_column(var = "gene_id") |>
  write_tsv(file.path(output_folder, "isoform_length_factors.tsv"))
