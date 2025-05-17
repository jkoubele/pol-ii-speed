library(GenomicFeatures)
library(AnnotationDbi)
library(txdbmaker)
library(tidyverse)


genome_folder <- "/home/jakub/Desktop/data_pol_ii/reference_genomes"
gtf_file_name <- "Homo_sapiens.GRCh38.112.gtf"

# Create TxDb object from your Ensembl GTF file
txdb <- makeTxDbFromGFF(file.path(genome_folder, gtf_file_name),
format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

write_tsv(tx2gene, file.path(genome_folder, 'tx2gene.tsv'))