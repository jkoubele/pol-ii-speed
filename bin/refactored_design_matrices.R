#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument(
  "--samplesheet",
  required = TRUE,
  help = "A .csv file with sample annotation. Must include a 'sample' column."
)
parser$add_argument(
  "--formula",
  required = TRUE,
  help = "Design formula using standard R syntax (e.g. '~ age + genotype')."
)


parser$add_argument(
  "--reduced_formulas",
  nargs = "+",
  help = "Optional list of reduced formulas, e.g. \"~ genotype\" \"~ age\"",
  default = NULL
)

parser$add_argument(
  "--pairwise_comparisons",
  nargs = "+",
  help = "NOT IMPLEMENTED YET. Intended format: factor,level1,level2 ...",
  default = NULL
)

parser$add_argument(
  "--output_folder",
  help = "",
  default = '.'
)

if (interactive()) {
  args <- list(
    samplesheet = "/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/design_matrix_test/samplesheet_test.csv",
    formula = "~ age + genotype",
    reduced_formulas = c("~ genotype", "~ age"),
    pairwise_comparisons = NULL,
    output_folder = "/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/design_matrix_test/output"
  )
} else {
  args <- parser$parse_args()
}

output_folder <- args$output_folder
reduced_matrices_output_folder <- file.path(output_folder, "reduced_design_matrices")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(reduced_matrices_output_folder, recursive = TRUE, showWarnings = FALSE)


samplesheet <- readr::read_csv(args$samplesheet)

## ---- Full design matrix ----
design_formula <- as.formula(args$formula)
terms_full <- terms(design_formula)
full_term_labels <- attr(terms_full, "term.labels")

X_full <- model.matrix(design_formula, samplesheet)
rownames(X_full) <- samplesheet$sample
X_full <- X_full[, colnames(X_full) != "(Intercept)", drop = FALSE]

df_full <- X_full |>
  as.data.frame() |>
  rownames_to_column(var = "sample")

write_csv(df_full, file.path(output_folder, "design_matrix_full.csv"))

## ---- Reduced formulas handling ----
# We'll collect metadata about each LRT here:
lrt_info <- tibble(
  name = character(),
  formula = character(),
  dropped_terms = character(),
  df_lrt = integer(),
  single_lfc = logical(),
  lfc_term = character()
)

# Build reduced_df from CLI (or empty if none provided)
if (is.null(args$reduced_formulas) || length(args$reduced_formulas) == 0) {
  reduced_df <- tibble()  # no LRTs specified
} else {
  reduced_df <- tibble(
    name = paste0("lrt_", seq_along(args$reduced_formulas)),
    formula = args$reduced_formulas
  )
}

if (nrow(reduced_df) > 0) {
  for (i in seq_len(nrow(reduced_df))) {
    test_name <- reduced_df$name[i]
    reduced_formula_string     <- reduced_df$formula[i]
    reduced_formula     <- as.formula(reduced_formula_string)

    # Check that reduced formula only drops terms (no new terms are added).
    terms_reduced <- terms(reduced_formula)
    terms_reduced_labels <- attr(terms_reduced, "term.labels")

    extra_terms <- setdiff(terms_reduced_labels, full_term_labels)
    if (length(extra_terms) > 0) {
      stop(sprintf(
        "Reduced formula '%s' contains terms not present in the full formula '%s': %s",
        reduced_formula_string, args$formula, paste(extra_terms, collapse = ", ")
      ))
    }

    dropped_terms <- setdiff(full_term_labels, terms_reduced_labels)

    if (length(dropped_terms) == 0) {
      stop(sprintf(
        "Reduced formula '%s' does not drop any term from the full formula '%s'.",
        reduced_formula_string, args$formula
      ))
    }

    # Build reduced design matrix
    reduced_design_matrix <- model.matrix(reduced_formula, samplesheet)
    rownames(reduced_design_matrix) <- samplesheet$sample
    reduced_design_matrix <- reduced_design_matrix[, colnames(X_red) != "(Intercept)", drop = FALSE]

    reduced_design_matrix_df <- reduced_design_matrix |>
      as.data.frame() |>
      rownames_to_column(var = "sample")

    out_path <- file.path(reduced_matrices_output_folder,
                          paste0("design_matrix_", test_name, ".csv"))
    write_csv(reduced_design_matrix_df, out_path)

    # 3) Compute df of the LRT and check if it's a single-LFC test
    df_lrt <- ncol(X_full) - ncol(reduced_design_matrix)

    # Single LFC if:
    #  - exactly one term was dropped in the formula
    #  - and that term corresponds to exactly 1 column in the design (df_lrt == 1)
    is_single_lfc <- (length(dropped_terms) == 1 && df_lrt == 1)
    lfc_term <- if (is_single_lfc) dropped_terms[1] else NA_character_

    lrt_info <- add_row(
      lrt_info,
      name = test_name,
      formula = reduced_formula_string,
      dropped_terms = paste(dropped_terms, collapse = ";"),
      df_lrt = df_lrt,
      single_lfc = is_single_lfc,
      lfc_term = lfc_term
    )
  }
}

write_csv(lrt_info, file.path(output_folder, "lrt_tests_metadata.csv"))

