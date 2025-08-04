#!/usr/bin/env Rscript

library(argparse)
library(rtracklayer)
library(dplyr)
library(readr)

parser <- ArgumentParser()
parser$add_argument("--gtf_file",
                    required = TRUE,
                    help = "Path to input GTF file.")
parser$add_argument("--output_folder",
                    default = '.')
args <- parser$parse_args()

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

gtf_df <- readGFF(args$gtf_file)

gene_type_column <- if ("gene_biotype" %in% names(gtf_df)) {
  "gene_biotype"
} else if ("gene_type" %in% names(gtf_df)) {
  "gene_type"
} else {
  stop("Neither 'gene_biotype' nor 'gene_type' column found in the GTF attributes.")
}

all_genes <- gtf_df |>
  dplyr::filter(type == "gene") |>
  dplyr::distinct()
protein_coding_genes <- all_genes |>
  dplyr::filter(!!sym(gene_type_column) == 'protein_coding') |>
  dplyr::select(gene_id)

readr::write_csv(all_genes |>  dplyr::select(gene_id), file.path(output_folder, "all_genes.csv"))
readr::write_csv(protein_coding_genes, file.path(output_folder, "protein_coding_genes.csv"))
