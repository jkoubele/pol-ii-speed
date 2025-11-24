#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(jsonlite)

parser <- ArgumentParser()
parser$add_argument("--samplesheet", required = TRUE,
                    help = "CSV with sample annotation. Must include a 'sample' column.")
parser$add_argument("--formula", required = TRUE,
                    help = "Design formula using standard R syntax (e.g. '~ age + genotype').")
parser$add_argument("--lrt_tests", default = NULL,
                    help = "Optional JSON file with LRT tests parsed from YAML.")
parser$add_argument("--output_folder", default = ".",
                    help = "Where to write design matrices and metadata.")

if (interactive()) {
  args <- list(
    samplesheet = "/home/jakub/Desktop/pol-ii-speed/design_matrix_test/samplesheet_test.csv",
    formula = "~ age + genotype + dummy_continuous",
    lrt_tests = "/home/jakub/Desktop/pol-ii-speed/design_matrix_test/lrt_tests.json",
    output_folder = "/home/jakub/Desktop/pol-ii-speed/design_matrix_test/output"
  )
} else {
  args <- parser$parse_args()
}

output_folder <- args$output_folder
output_folder_reduced_matrices <- file.path(output_folder, "reduced_design_matrices")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(output_folder_reduced_matrices, recursive = TRUE, showWarnings = FALSE)

samplesheet <- read_csv(args$samplesheet, show_col_types = FALSE)
if (!"sample" %in% names(samplesheet)) {
  stop("samplesheet must contain a column named 'sample'.")
}

lrt_tests <- list()
if (!is.null(args$lrt_tests)) {
  lrt_tests <- fromJSON(args$lrt_tests, simplifyVector = FALSE)
}

design_formula <- as.formula(args$formula)
terms_full <- terms(design_formula)
full_term_labels <- attr(terms_full, "term.labels")

design_matrix <- model.matrix(design_formula, samplesheet)
rownames(design_matrix) <- samplesheet$sample
design_matrix <- design_matrix[, colnames(design_matrix) != "(Intercept)", drop = FALSE]

design_matrix_df <- design_matrix |>
  as.data.frame() |>
  rownames_to_column("sample")

write_csv(design_matrix_df, file.path(output_folder, "design_matrix.csv"))

lrt_metadata <- tibble(
  test_name = character(),
  variable = character(),
  test_type = character(),  # "continuous" or "categorical"
  group_1 = character(),
  group_2 = character(),
  lrt_df = integer(),
  lfc_column_positive = character(),  # optional
  lfc_column_negative = character()   # optional
)

if (length(lrt_tests) == 0) {
  write_csv(lrt_metadata, file.path(output_folder, "lrt_tests_metadata.csv"))
  warning("No LRT specification provided, creating only design matrix to estimate effect sizes.")
  quit(status = 0)
}

for (i in seq_along(lrt_tests)) {
  test_dict <- lrt_tests[[i]]
  variable_name <- as.character(test_dict$variable)
  group_1 <- if ("group_1" %in% names(test_dict)) as.character(test_dict$group_1) else NULL
  group_2 <- if ("group_2" %in% names(test_dict)) as.character(test_dict$group_2) else NULL

  if (!variable_name %in% names(samplesheet)) {
    stop(sprintf("Variable '%s' not found in samplesheet columns.", variable_name))
  }

  test_name <- paste0("lrt_", i)
  reduced_matrix_output_path <- file.path(
    output_folder_reduced_matrices,
    paste0(test_name, ".csv")
  )

  lfc_column_positive <- ""
  lfc_column_negative <- ""

  if (is.numeric(samplesheet[[variable_name]])) {
    test_type <- "continuous"
    
    if (!is.null(group_1) || !is.null(group_2)) {
      stop(sprintf("Continuous variable '%s' cannot specify comparison groups.", variable_name))
    }
    stopifnot(variable_name %in% full_term_labels)

    reduced_terms <- setdiff(full_term_labels, variable_name)
    reduced_formula_str <- if (length(reduced_terms) == 0) "~ 1" else paste("~", paste(reduced_terms, collapse = " + "))
    reduced_formula <- as.formula(reduced_formula_str)

    reduced_design_matrix <- model.matrix(reduced_formula, samplesheet)
    rownames(reduced_design_matrix) <- samplesheet$sample
    reduced_design_matrix <- reduced_design_matrix[, colnames(reduced_design_matrix) != "(Intercept)", drop = FALSE]

    lrt_df <- ncol(design_matrix) - ncol(reduced_design_matrix)
    lfc_column_positive <- variable_name
    lfc_column_negative <- ""

    reduced_design_matrix_df <- reduced_design_matrix |>
      as.data.frame() |>
      rownames_to_column("sample")
    write_csv(reduced_design_matrix_df, reduced_matrix_output_path)

    group_1 <- ""
    group_2 <- ""

  } else {
    # Categorical test
    test_type <- "categorical"
    stopifnot(variable_name %in% full_term_labels)

    samplesheet[[variable_name]] <- factor(samplesheet[[variable_name]])
    levels_full <- levels(samplesheet[[variable_name]])
    num_levels <- length(levels_full)
    stopifnot(num_levels >= 2)

    # Auto-generate groups if both are missing (for binary variables only)
    if (is.null(group_1) && is.null(group_2)) {
      if (num_levels != 2) {
        stop(sprintf(
          "Categorical variable '%s' has %d levels. You must specify group_1 and group_2.",
          variable_name, num_levels
        ))
      }
      group_1 <- levels_full[1]
      group_2 <- levels_full[2]
    } else if (is.null(group_1) || is.null(group_2)) {
      stop("Only one comparison group was specified.")
    }

    stopifnot(group_1 != group_2)
    stopifnot(group_1 %in% levels_full)
    stopifnot(group_2 %in% levels_full)

    baseline_full <- levels_full[1]

    dummy_matrix_full <- model.matrix(reformulate(variable_name), samplesheet)
    dummy_matrix_full <- dummy_matrix_full[, colnames(dummy_matrix_full) != "(Intercept)", drop = FALSE]
    nonbase_levels_full <- levels_full[-1]
    level_to_column_full <- setNames(colnames(dummy_matrix_full), nonbase_levels_full)

    if (group_1 != baseline_full) {
      lfc_column_positive <- level_to_column_full[[group_1]]
      if (lfc_column_positive != "") {
        stopifnot(!is.na(lfc_column_positive), lfc_column_positive %in% colnames(design_matrix))
      }
    } else {
      lfc_column_positive <- ""
    }

    if (group_2 != baseline_full) {
      lfc_column_negative <- level_to_column_full[[group_2]]
      if (lfc_column_negative != "") {
        stopifnot(!is.na(lfc_column_negative), lfc_column_negative %in% colnames(design_matrix))
      }
    } else {
      lfc_column_negative <- ""
    }

    # Build reduced matrix with df=1 by releveling to group_2 baseline
    samplesheet_relevelled <- samplesheet
    samplesheet_relevelled[[variable_name]] <- relevel(samplesheet_relevelled[[variable_name]], ref = group_2)

    full_matrix_relevelled <- model.matrix(design_formula, samplesheet_relevelled)
    rownames(full_matrix_relevelled) <- samplesheet_relevelled$sample
    full_matrix_relevelled <- full_matrix_relevelled[, colnames(full_matrix_relevelled) != "(Intercept)", drop = FALSE]

    dummy_matrix_relevelled <- model.matrix(reformulate(variable_name), samplesheet_relevelled)
    dummy_matrix_relevelled <- dummy_matrix_relevelled[, colnames(dummy_matrix_relevelled) != "(Intercept)", drop = FALSE]
    levels_relevelled <- levels(samplesheet_relevelled[[variable_name]])
    nonbase_relevelled <- levels_relevelled[-1]
    level_to_col_relevelled <- setNames(colnames(dummy_matrix_relevelled), nonbase_relevelled)

    col_drop <- level_to_col_relevelled[[group_1]]
    if (is.na(col_drop) || col_drop == "") {
      stop(sprintf("Could not find contrast column for %s: %s vs %s.",
                   variable_name, group_1, group_2))
    }
    stopifnot(col_drop %in% colnames(full_matrix_relevelled))
    
    reduced_design_matrix <- full_matrix_relevelled[, colnames(full_matrix_relevelled) != col_drop, drop = FALSE]
    lrt_df <- ncol(full_matrix_relevelled) - ncol(reduced_design_matrix)
    stopifnot(lrt_df == 1)

    reduced_design_matrix_df <- reduced_design_matrix |>
      as.data.frame() |>
      rownames_to_column("sample")
    write_csv(reduced_design_matrix_df, reduced_matrix_output_path)
  }

  lrt_metadata <- add_row(
    lrt_metadata,
    test_name = test_name,
    variable = variable_name,
    test_type = test_type,
    group_1 = group_1,
    group_2 = group_2,
    lrt_df = lrt_df,
    lfc_column_positive = lfc_column_positive,
    lfc_column_negative = lfc_column_negative
  )
  
}

write_csv(lrt_metadata, file.path(output_folder, "lrt_tests_metadata.csv"))
