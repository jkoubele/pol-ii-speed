#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)


parser <- ArgumentParser()
parser$add_argument("--regularization_chunks_folder",
                    required = TRUE,
                    help = "Folder containing CSV files with regularization chunk results.")
parser$add_argument("--test_results",
                    required = TRUE,
                    help = "CSV file with LRT results.")
parser$add_argument("--output_folder",
                    default = '.',
                    help = 'Path to output folder')

args <- parser$parse_args()

output_folder <- args$output_folder

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

regularized_model_parameters_files <- list.files(path = args$regularization_chunks_folder,
                                     pattern = "^regularized_model_parameters.*\\.csv$",
                                     full.names = TRUE)

if (length(regularized_model_parameters_files) == 0) {
  stop("No matching CSV files with regularized model parameters found in the folder: ", args$regularization_chunks_folder)
}

regularized_model_parameters_merged_df <- sort(regularized_model_parameters_files) |>
  map(read_csv) |>
  keep(~nrow(.x) > 0) |>
  list_rbind()

write_csv(regularized_model_parameters_merged_df, file.path(output_folder, 'regularized_model_parameters.csv'))

test_results <- read_csv(args$test_results)

lfc_value_positive <- test_results |>
  left_join(
    regularized_model_parameters_merged_df,
    by = c(
      "tested_parameter"   = "parameter_type",
      "gene_name"          = "gene_name",
      "intron_name"        = "intron_name",
      "lfc_column_positive" = "feature_name"
    )
  ) |>
  mutate(value = replace_na(value, 0.0)) |>
  pull(value)

lfc_value_negative <- test_results |>
  left_join(
    regularized_model_parameters_merged_df,
    by = c(
      "tested_parameter"   = "parameter_type",
      "gene_name"          = "gene_name",
      "intron_name"        = "intron_name",
      "lfc_column_negative" = "feature_name"
    )
  ) |>
  mutate(value = replace_na(value, 0.0)) |>
  pull(value)

test_results$lfc_regularized <- lfc_value_positive - lfc_value_negative

write_csv(test_results, file.path(output_folder, 'test_results.csv'))
write_csv(test_results, file.path(output_folder, 'test_results_full.csv'))
