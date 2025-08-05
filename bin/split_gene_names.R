#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)


parser <- ArgumentParser()
parser$add_argument("--input_gene_names", required = TRUE, help = "Path to input CSV file with gene names.")
parser$add_argument("--output_folder", required = TRUE, help = "Directory to write output chunked CSV files.")
parser$add_argument("--chunk_size", type = "integer", default = 100, help = "Number of rows per chunk (default: 100).")
args <- parser$parse_args()

input_file <- args$input_gene_names
output_folder <- args$output_folder
chunk_size <- args$chunk_size

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

gene_names <- read_csv(input_file)

n_chunks <- ceiling(nrow(gene_names) / chunk_size)
pad_width <- nchar(as.character(n_chunks))  # Determine left-padding


for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, nrow(gene_names))

  chunk <- gene_names[start:end,]

  filename <- sprintf("gene_names_chunk_%0*d.csv", pad_width, i)
  output_path <- file.path(output_folder, filename)

  write_csv(chunk, output_path)
}

