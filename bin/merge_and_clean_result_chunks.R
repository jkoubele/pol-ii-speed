#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)


parser <- ArgumentParser()
parser$add_argument("--input_folder",
                    required = TRUE,
                    help = "Folder containing CSV files to merge.")
parser$add_argument("--output_folder",
                    default = '.',
                    help = 'Path to output folder')

args <- parser$parse_args()

input_folder <- args$input_folder
output_folder <- args$output_folder

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

model_parameters_files <- list.files(path = input_folder,
                                     pattern = "^model_parameters.*\\.csv$",
                                     full.names = TRUE)

if (length(model_parameters_files) == 0) {
  stop("No matching CSV files with model parameters found in the input folder: ", input_folder)
}

model_parameters_merged_df <- sort(model_parameters_files) |>
  map(read_csv) |>
  keep(~nrow(.x) > 0) |>
  list_rbind()

write_csv(model_parameters_merged_df, file.path(output_folder, 'model_parameters.csv'))


test_results_files <- list.files(path = input_folder,
                                 pattern = "^test_results.*\\.csv$",
                                 full.names = TRUE)

if (length(test_results_files) == 0) {
  stop("No matching CSV files with test results found in the input folder: ", input_folder)
}

test_results_merged_df <- sort(test_results_files) |>
  map(read_csv) |>
  keep(~nrow(.x) > 0) |>
  list_rbind()

write_csv(test_results_merged_df, file.path(output_folder, 'test_results_raw.csv'))

test_results_cleaned_df <- test_results_merged_df |>
  mutate(
    parameter_type = case_when(
      tested_parameter == "alpha" ~ "gene_expression",
      tested_parameter == "beta" ~ "elongation_speed",
      tested_parameter == "gamma" ~ "splicing_speed",
      TRUE ~ tested_parameter
    ),
    # replace splicing time (parameter gamma) by splicing speed
    lfc = if_else(parameter_type == "splicing_speed", -lfc, lfc),
    l2fc = lfc / log(2)
  ) |>
  select(gene_name, intron_name, parameter_type, test_id, variable, group_1, group_2, l2fc, p_value, chi2_test_statistics)

write_csv(test_results_cleaned_df, file.path(output_folder, 'test_results.csv'))


ignored_introns_files <- sort(list.files(path = input_folder,
                                         pattern = "^ignored_introns.*\\.csv$",
                                         full.names = TRUE))

if (length(ignored_introns_files) > 0) {
  ignored_introns_files |>
    map(~read_csv(.x, show_col_types = FALSE)) |>
    list_rbind() |>
    write_csv(file.path(output_folder, "ignored_introns.csv"))
} else {
  message("No ignored_introns*.csv files found in input folder: ", input_folder)
}

ignored_genes_files <- sort(list.files(path = input_folder,
                                       pattern = "^ignored_genes.*\\.csv$",
                                       full.names = TRUE))

if (length(ignored_genes_files) > 0) {
  ignored_genes_files |>
    map(~read_csv(.x, show_col_types = FALSE)) |>
    list_rbind() |>
    write_csv(file.path(output_folder, "ignored_genes.csv"))
} else {
  message("No ignored_genes*.csv files found in input folder: ", input_folder)
}
