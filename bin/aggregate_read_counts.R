#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(tximport)
library(DESeq2)


parser <- ArgumentParser()
parser$add_argument("--tx2gene",
                    required = TRUE,
                    help = "Path to tx2gene.tsv")
parser$add_argument("--exon_quant_files", nargs = "+", required = TRUE)
parser$add_argument("--intron_counts_files", nargs = "+", required = TRUE)
parser$add_argument("--sample_names", nargs = "+", required = TRUE)
parser$add_argument("--output_folder", default = '.')

args <- parser$parse_args()

output_folder <- args$output_folder
sample_names <- args$sample_names
exon_quant_files <- args$exon_quant_files
intron_counts_files <- args$intron_counts_files

tx2gene <- read_tsv(args$tx2gene, show_col_types = FALSE) |>
  mutate(
    TXNAME = sub("\\|.*$", "", TXNAME),
    TXNAME = sub("\\.[0-9]+$", "", TXNAME)
)



if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

stopifnot(length(sample_names) == length(exon_quant_files),
          length(sample_names) == length(intron_counts_files))


names(exon_quant_files) <- args$sample_names
txi <- tximport(exon_quant_files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE,
                countsFromAbundance = 'scaledTPM')

exon_read_counts <- as.data.frame(txi$counts) |>
  rownames_to_column(var = "gene_id")

intron_counts_by_sample <- map2(sample_names, intron_counts_files, function(sample, file) {
  read_tsv(file) |>
    select(name, count) |>
    set_names(c("name", sample))
})

intron_read_counts <- purrr::reduce(intron_counts_by_sample, full_join, by = "name") |>
  dplyr::rename(intron_id = name)


exon_read_counts <- exon_read_counts |> select(gene_id, all_of(sample_names))
intron_read_counts <- intron_read_counts |> select(intron_id, all_of(sample_names))

write_tsv(exon_read_counts, file.path(args$output_folder, 'exon_counts.tsv'))
write_tsv(intron_read_counts, file.path(args$output_folder, 'intron_counts.tsv'))

# Estimate library sizes by DESeq2 method
all_reads_matrix <- bind_rows(column_to_rownames(exon_read_counts, var = 'gene_id'),
                              column_to_rownames(intron_read_counts, var = 'intron_id')) |> round()

dds <- DESeqDataSetFromMatrix(countData = all_reads_matrix,
                              colData = DataFrame(row.names = colnames(all_reads_matrix)),
                              design = ~1) |>
  estimateSizeFactors()

size_factors <- sizeFactors(dds) |>
  enframe(name = "sample", value = "library_size_factor")
write_tsv(size_factors, file.path(output_folder, "library_size_factors.tsv"))
