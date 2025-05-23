library(argparse)
library(tidyverse)
library(tximport)
library(DESeq2)


parser <- ArgumentParser()
parser$add_argument("--tx2gene",
                    help = "Path to tx2gene.tsv",
                    default = '/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/tx2gene.tsv')
parser$add_argument("--salmon_folder",
                    help = "",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/salmon_output')
parser$add_argument("--intronic_reads_folder",
                    help = "",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intronic_reads')
parser$add_argument("--output_folder",
                    help = "",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/aggregated_counts')
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
                ignoreTxVersion = FALSE,
                countsFromAbundance = 'scaledTPM')

exon_reads_counts <- as.data.frame(txi$counts) |>
  rownames_to_column(var = "gene_id")

dir.create(args$output_folder, recursive = TRUE, showWarnings = FALSE)
write_tsv(exon_reads_counts, file.path(args$output_folder, 'exonic_reads.tsv'))


intron_counts_by_sample <- list()
for (sample_name in list.dirs(args$intronic_reads_folder, recursive = FALSE, full.names = FALSE)) {
  df_intron_counts <- read_tsv(file.path(args$intronic_reads_folder, sample_name, "intron_read_counts.tsv")) |>
    select(name, count) |>
    rename_with(~sample_name, .cols = "count")
  intron_counts_by_sample[[sample_name]] <- df_intron_counts
}

intron_counts <- purrr::reduce(intron_counts_by_sample, full_join, by = "name")
write_tsv(intron_counts, file.path(args$output_folder, 'intron_counts.tsv'))

# Estimate library sizes by DESeq2 method
all_reads_matrix <- bind_rows(column_to_rownames(exon_reads_counts, var = 'gene_id'),
                              column_to_rownames(intron_counts, var = 'name')) |>
  round()

dds <- DESeqDataSetFromMatrix(countData = all_reads_matrix,
                              colData = DataFrame(row.names = colnames(all_reads_matrix)),
                              design = ~1) |>
  estimateSizeFactors()

sizeFactors(dds) |>
  enframe(name = "sample_name", value = "library_size_factor") |>
  write_tsv(file.path(args$output_folder, "library_size_factors.tsv"))

