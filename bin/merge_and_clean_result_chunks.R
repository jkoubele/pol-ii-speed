#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)


parser <- ArgumentParser()
parser$add_argument("--input_folder",
                    required = TRUE,
                    help = "Folder containing CSV files to merge.")
parser$add_argument("--output_folder", 
                    default = '.',
                    help='Path to output folder')

args <- parser$parse_args()

input_folder <- args$input_folder
output_folder <- args$output_folder

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

input_files <- list.files(path = args$input_folder,
                          pattern = "\\.csv$",
                          full.names = TRUE)

if (length(input_files) == 0) {
  stop("No matching CSV files found in input folder: ", args$input_folder)
}

merged_df <- sort(input_files) |>
  map(read_csv) |>
  keep(~ nrow(.x) > 0) |>
  list_rbind()

write_csv(merged_df, file.path(output_folder, 'model_results_raw.csv'))

cleaned_df <- merged_df |>
  mutate(
    chi2_test_statistic = 2 * loss_differences,
    parameter_type = case_when(
      parameter_type == "alpha" ~ "gene_expression",
      parameter_type == "beta"  ~ "elongation_speed",
      parameter_type == "gamma" ~ "splicing_speed",
      TRUE ~ parameter_type
    ),
    # replace splicing time (parameter gamma) by splicing speed
    value = if_else(parameter_type == "splicing_speed", -value, value),
    log2_fold_change = value / log(2),
    p_value = p_value_lrt
  ) |>
  filter(parameter_type %in% c("gene_expression", "elongation_speed", "splicing_speed")) |>
  select(parameter_type, log2_fold_change, feature_name, gene_name, intron_name, chi2_test_statistic, p_value)

if (all(is.na(cleaned_df$intron_name))) {
  cleaned_df <- cleaned_df |> select(-intron_name)
}

write_csv(cleaned_df, file.path(output_folder, 'model_results.csv'))
