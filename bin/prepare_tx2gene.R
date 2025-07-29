#!/usr/bin/env Rscript

library(argparse)
library(GenomicFeatures)
library(AnnotationDbi)
library(readr)


parser <- ArgumentParser()
parser$add_argument("--gtf_file",
                    required = TRUE,
                    help = "Path to input GTF file.")
parser$add_argument("--output_folder",
                    default = '.')
args <- parser$parse_args()


if (!dir.exists(args$output_folder)) {
  dir.create(args$output_folder, recursive = TRUE)
}

txdb <- makeTxDbFromGFF(file.path(args$gtf_file), format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

write_tsv(tx2gene, file.path(args$output_folder, 'tx2gene.tsv'))