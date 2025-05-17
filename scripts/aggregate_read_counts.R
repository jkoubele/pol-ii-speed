library(argparse)
library(tidyverse)
library(tximport)


parser <- ArgumentParser()
parser$add_argument("--tx2gene",
                    help = "Path to tx2gene.tsv",
                    default = '/home/jakub/Desktop/data_pol_ii/reference_genomes/tx2gene.tsv')
parser$add_argument("--salmon_folder",
                    help = "",
                    default = '/home/jakub/Desktop/data_pol_ii/human_astrocytes/salmon_output')
parser$add_argument("--intronic_reads_folder",
                    help = "",
                    default = '/home/jakub/Desktop/data_pol_ii/human_astrocytes/intronic_reads')
parser$add_argument("--output_folder",
                    help = "",
                    default = '/home/jakub/Desktop/data_pol_ii/human_astrocytes/aggregated_counts')
args <- parser$parse_args()

tx2gene <- read_tsv(args$tx2gene)

files <- c()
for (sample_name in list.dirs(args$salmon_folder, recursive = FALSE, full.names = FALSE)) {
  file_path <- file.path(args$salmon_folder, sample_name, "quant.sf")
  files[sample_name] <- file_path
}

txi <- tximport(files,
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion=TRUE,
                countsFromAbundance='scaledTPM')

exon_reads_counts <- as.data.frame(txi$counts) |>
  rownames_to_column(var = "gene_id")

dir.create(args$output_folder, recursive = TRUE, showWarnings = FALSE)
write_tsv(exon_reads_counts, file.path(args$output_folder, 'exonic_reads.tsv'))


intron_counts_by_sample <- list()
for (sample_name in list.dirs(args$intronic_reads_folder, recursive = FALSE, full.names = FALSE)) {
  df_intron_counts <- read_tsv(file.path(args$intronic_reads_folder, sample_name, "intron_read_counts.tsv")) |>
    select(name, count) |>
    rename(!!sample_name := count)
  intron_counts_by_sample[[sample_name]] <- df_intron_counts
}

intron_counts <- reduce(intron_counts_by_sample, full_join, by = "name")
write_tsv(intron_counts, file.path(args$output_folder, 'intron_counts.tsv'))
