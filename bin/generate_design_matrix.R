#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("--samplesheet",
                    required = TRUE,
                    help = "A .csv file with sample annotation. Must include a 'sample' column.")
parser$add_argument("--formula",
                    required = TRUE,
                    help = "Design formula using standard R syntax (e.g. '~ age + genotype' or '~ age:condition')")
parser$add_argument("--factor_references",
                    help = "An optional .csv file with factor base levels.")
parser$add_argument("--output_folder",
                    help = "",
                    default = '.')

args <- parser$parse_args()

output_folder <- args$output_folder

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

samplesheet <- readr::read_csv(args$samplesheet)

if (!is.null(args$factor_references)) {
  factor_references <- readr::read_csv(args$factor_references)
  for (i in seq_len(nrow(factor_references))) {
    var <- factor_references$factor_name[i]
    ref <- factor_references$reference[i]

    if (var %in% colnames(samplesheet)) {
      samplesheet[[var]] <- forcats::fct_relevel(as.factor(samplesheet[[var]]), ref)
    } else {
      warning(sprintf("Column '%s' specified in factor_references not found in samplesheet.", var))
    }
  }

}

design_formula <- as.formula(args$formula)
design_matrix <- model.matrix(design_formula, samplesheet)

rownames(design_matrix) <- samplesheet$sample
design_matrix <- design_matrix[, colnames(design_matrix) != "(Intercept)", drop = FALSE]

df_out <- design_matrix |>
  as.data.frame() |>
  rownames_to_column(var = "sample")

write_csv(df_out, file.path(output_folder, 'design_matrix.csv'))
