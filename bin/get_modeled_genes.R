#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)


parser <- ArgumentParser()
parser$add_argument("--input_gene_names", required = TRUE, help = "Path to input CSV file with gene names.")
parser$add_argument("--modelable_genes", required = TRUE, help = "Path to TSV file with modelable genes.")
parser$add_argument("--modelable_introns", required = TRUE, help = "Path to TSV file with modelable introns.")
parser$add_argument("--output_folder", required = TRUE, help = "Directory to write output chunked CSV files.")
parser$add_argument("--chunk_size", type = "integer", default = 100, help = "Number of rows per chunk (default: 100).")
args <- parser$parse_args()


output_folder <- args$output_folder
chunk_size <- args$chunk_size

input_gene_names <- read_csv(args$input_gene_names)
modelable_genes <- read_tsv(args$modelable_genes)
modelable_introns <- read_tsv(args$modelable_introns)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

chunks_output_subfolder = file.path(output_folder, '/gene_chunks')
if (!dir.exists(chunks_output_subfolder)) {
  dir.create(chunks_output_subfolder, recursive = TRUE)
}

modeled_genes <- modelable_genes |>
  semi_join(input_gene_names, by = "gene_id")

modeled_introns <- modelable_introns |>
  semi_join(input_gene_names, by = "gene_id")

modeled_genes |>
  write_tsv(file.path(output_folder, "modeled_genes.tsv"))

modeled_introns |>
  write_tsv(file.path(output_folder, "modeled_introns.tsv"))


n_chunks <- ceiling(nrow(modeled_genes) / chunk_size)
pad_width <- nchar(as.character(n_chunks))  # Determine left-padding

for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, nrow(modeled_genes))

  chunk <- modeled_genes[start:end,]

  filename <- sprintf("modeled_genes_chunk_%0*d.tsv", pad_width, i)
  output_path <- file.path(chunks_output_subfolder, filename)

  write_tsv(chunk, output_path)
}
