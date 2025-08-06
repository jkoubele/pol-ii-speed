#!/usr/bin/env Rscript

library(argparse)
library(readr)
library(purrr)


parser <- ArgumentParser()
parser$add_argument("--input_folder",
                    required = TRUE,
                    help = "Folder containing CSV files to merge.")
parser$add_argument("--output_folder", 
                    default = '.',
                    help='Path to output folder')
parser$add_argument("--output_file_name", 
                    default = 'model_results.csv',
                    help = "Name of output CSV file.")

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
  map_dfr(read_csv)

write_csv(merged_df, file.path(output_folder, args$output_file_name))
