library(argparse)
library(GenomicFeatures)
library(AnnotationDbi)
library(txdbmaker)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("--gtf_file",
                    help = "",
                    default = '/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/Caenorhabditis_elegans.WBcel235.112.gtf')
parser$add_argument("--output_folder",
                    help = "",
                    default = '/cellfile/datapublic/jkoubele/reference_genomes/WBcel235')
args <- parser$parse_args()


# Create TxDb object from your Ensembl GTF file
txdb <- makeTxDbFromGFF(file.path(args$gtf_file), format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

write_tsv(tx2gene, file.path(args$output_folder, 'tx2gene.tsv'))