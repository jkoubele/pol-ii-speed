#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(jsonlite)


parser <- ArgumentParser()
parser$add_argument("--samplesheet", required=TRUE,
                    help="CSV with sample annotation. Must include a 'sample' column.")
parser$add_argument("--formula", required=TRUE,
                    help="Design formula using standard R syntax (e.g. '~ age + genotype').")
parser$add_argument("--lrt_tests", default=NULL,
                    help="Optional JSON file with LRT tests parsed from YAML.")
parser$add_argument("--output_folder", default=".",
                    help="Where to write design matrices and metadata.")

if (interactive()) {
  args <- list(
    samplesheet   = "/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/design_matrix_test/samplesheet_test.csv",
    formula       = "~ age + genotype",
    lrt_tests     = "/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/design_matrix_test/lrt_tests.json",
    output_folder = "/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/design_matrix_test/output"
  )
} else {
  args <- parser$parse_args()
}

outdir <- args$output_folder
reduced_dir <- file.path(outdir, "reduced_design_matrices")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(reduced_dir, recursive=TRUE, showWarnings=FALSE)

# -------------- Read samplesheet --------------

samplesheet <- readr::read_csv(args$samplesheet, show_col_types = FALSE)

if (!"sample" %in% names(samplesheet)) {
  stop("samplesheet must contain a column named 'sample'.")
}

# -------------- Parse tests JSON --------------

lrt_tests <- list()
if (!is.null(args$lrt_tests)) {
  if (!file.exists(args$lrt_tests)) {
    stop(sprintf("lrt_tests JSON file '%s' not found.", args$lrt_tests))
  }
  lrt_tests <- jsonlite::fromJSON(args$lrt_tests, simplifyVector = FALSE)
}

# -------------- Helpers --------------

# Safe getter for optional fields in JSON objects
get_or_empty <- function(x) {
  if (is.null(x) || is.na(x) || x == "") "" else as.character(x)
}

# Extract columns in a model.matrix corresponding to a single *standalone* term
# using the "assign" attribute.
term_columns <- function(X, formula_obj, term_label) {
  tt <- terms(formula_obj)
  labels <- attr(tt, "term.labels")
  
  idx <- match(term_label, labels)
  if (is.na(idx)) return(character())
  
  assign_vec <- attr(X, "assign")  # 0 is intercept, 1..k correspond to term.labels
  cols <- colnames(X)[assign_vec == idx]
  cols
}

# Ensure deterministic factor baselines for the FULL matrix:
# for every plain variable used in the formula that is categorical,
# reorder levels alphabetically.
design_formula <- as.formula(args$formula)
terms_full <- terms(design_formula)
full_term_labels <- attr(terms_full, "term.labels")

for (term in full_term_labels) {
  # only plain variable names, not interactions or transforms
  if (grepl("[:*()^]|I\\(", term)) next
  if (!term %in% names(samplesheet)) next
  
  if (is.character(samplesheet[[term]]) || is.factor(samplesheet[[term]])) {
    levs <- sort(unique(as.character(samplesheet[[term]])))
    samplesheet[[term]] <- factor(samplesheet[[term]], levels = levs)
  }
}

# -------------- Build FULL matrix once --------------

X_full <- model.matrix(design_formula, samplesheet)
rownames(X_full) <- samplesheet$sample
X_full <- X_full[, colnames(X_full) != "(Intercept)", drop = FALSE]

df_full <- X_full |>
  as.data.frame() |>
  rownames_to_column("sample")

# Backwards-compatible name expected by your pipeline
readr::write_csv(df_full, file.path(outdir, "design_matrix.csv"))
# Optional clarity copy
readr::write_csv(df_full, file.path(outdir, "design_matrix_full.csv"))

# -------------- Prepare metadata tibble --------------

lrt_info <- tibble(
  test_name        = character(),
  variable         = character(),
  test_type        = character(),  # "continuous" or "categorical"
  group1           = character(),
  group2           = character(),
  restricted_matrix= character(),
  df_lrt           = integer(),
  lfc_col1         = character(),  # optional
  lfc_col2         = character()   # optional
)

# If no tests, still write empty metadata
if (length(lrt_tests) == 0) {
  readr::write_csv(lrt_info, file.path(outdir, "lrt_tests_metadata.csv"))
  quit(status = 0)
}

# -------------- Main loop over tests --------------

for (i in seq_along(lrt_tests)) {
  
  test <- lrt_tests[[i]]
  var <- get_or_empty(test$variable)
  g1  <- get_or_empty(test$group1)
  g2  <- get_or_empty(test$group2)
  
  if (var == "") stop("Each LRT test must include 'variable'.")
  
  if (!var %in% names(samplesheet)) {
    stop(sprintf("LRT variable '%s' not found in samplesheet columns.", var))
  }
  
  test_name <- paste0("lrt_", i)
  restricted_path <- file.path(reduced_dir, paste0("design_matrix_", test_name, ".csv"))
  
  xcol1 <- ""
  xcol2 <- ""
  
  # Decide continuous vs categorical
  is_numeric <- is.numeric(samplesheet[[var]])
  
  if (is_numeric) {
    # ---------- CONTINUOUS TEST ----------
    if (g1 != "" || g2 != "") {
      stop(sprintf("Continuous variable '%s' must not specify group1/group2.", var))
    }
    
    test_type <- "continuous"
    
    # Reduced formula = full formula without this standalone term
    if (!var %in% full_term_labels) {
      stop(sprintf(
        "Continuous test for '%s' but term is not a standalone term in full formula '%s'.",
        var, args$formula
      ))
    }
    
    reduced_terms <- setdiff(full_term_labels, var)
    reduced_formula_str <- if (length(reduced_terms) == 0) "~ 1"
    else paste("~", paste(reduced_terms, collapse=" + "))
    reduced_formula <- as.formula(reduced_formula_str)
    
    X_red <- model.matrix(reduced_formula, samplesheet)
    rownames(X_red) <- samplesheet$sample
    X_red <- X_red[, colnames(X_red) != "(Intercept)", drop = FALSE]
    
    df_lrt <- ncol(X_full) - ncol(X_red)
    if (df_lrt != 1) {
      stop(sprintf(
        "Continuous test '%s' produced df_lrt=%d (expected 1). Check formula/coding.",
        test_name, df_lrt
      ))
    }
    
    # LFC column(s) in FULL matrix for this term
    cols_full <- term_columns(X_full, design_formula, var)
    if (length(cols_full) != 1) {
      stop(sprintf(
        "Continuous term '%s' maps to %d columns in full matrix (expected 1).",
        var, length(cols_full)
      ))
    }
    xcol1 <- cols_full[1]
    xcol2 <- ""
    
    # Write reduced CSV
    df_red <- X_red |>
      as.data.frame() |>
      rownames_to_column("sample")
    readr::write_csv(df_red, restricted_path)
    
    group1_res <- ""
    group2_res <- ""
    
  } else {
    # ---------- CATEGORICAL TEST ----------
    test_type <- "categorical"
    
    samplesheet[[var]] <- factor(samplesheet[[var]])
    levs_full <- levels(samplesheet[[var]])
    nlev <- length(levs_full)
    if (nlev < 2) stop(sprintf("Variable '%s' has <2 levels.", var))
    
    # Auto-generate groups for binary if missing
    if (g1 == "" && g2 == "") {
      if (nlev != 2) {
        stop(sprintf(
          "Categorical variable '%s' has %d levels. You must specify group1 and group2 explicitly.",
          var, nlev
        ))
      }
      # deterministic direction: alphabetical levels already set for FULL matrix
      g2 <- levs_full[1]
      g1 <- levs_full[2]
    } else if (g1 == "" || g2 == "") {
      stop(sprintf(
        "Categorical test for '%s' must specify both group1 and group2, or neither (binary only).",
        var
      ))
    }
    
    if (g1 == g2) stop(sprintf("group1 and group2 must differ for '%s'.", var))
    if (!g1 %in% levs_full) stop(sprintf("group1 '%s' not a level of '%s'.", g1, var))
    if (!g2 %in% levs_full) stop(sprintf("group2 '%s' not a level of '%s'.", g2, var))
    
    # ---- LFC columns relative to FULL baseline ----
    baseline_full <- levs_full[1]
    
    # Build one-variable matrix to get column name mapping for FULL baseline
    mm_full <- model.matrix(reformulate(var), samplesheet)
    mm_full <- mm_full[, colnames(mm_full) != "(Intercept)", drop=FALSE]
    nonbase_levels_full <- levs_full[-1]
    level_to_col_full <- setNames(colnames(mm_full), nonbase_levels_full)
    
    # Determine lfc_col1 (group1 - baseline_full), optional if group1 is baseline
    if (g1 != baseline_full) {
      xcol1 <- level_to_col_full[[g1]]
      if (is.na(xcol1) || xcol1 == "") stop("Internal error mapping group1 to column.")
      if (!xcol1 %in% colnames(X_full)) stop("lfc_col1 not found in full matrix.")
    } else {
      xcol1 <- ""
    }
    
    # Determine lfc_col2 (group2 - baseline_full), optional if group2 is baseline
    if (g2 != baseline_full) {
      xcol2 <- level_to_col_full[[g2]]
      if (is.na(xcol2) || xcol2 == "") stop("Internal error mapping group2 to column.")
      if (!xcol2 %in% colnames(X_full)) stop("lfc_col2 not found in full matrix.")
    } else {
      xcol2 <- ""
    }
    
    # ---- Reduced matrix construction with df=1 ----
    # Relevel temporarily so group2 is baseline, then drop group1 column.
    
    ss2 <- samplesheet
    ss2[[var]] <- relevel(ss2[[var]], ref = g2)
    
    X_full2 <- model.matrix(design_formula, ss2)
    rownames(X_full2) <- ss2$sample
    X_full2 <- X_full2[, colnames(X_full2) != "(Intercept)", drop=FALSE]
    
    # Identify the single column encoding (group1 - group2) in this relevelled coding
    mm2 <- model.matrix(reformulate(var), ss2)
    mm2 <- mm2[, colnames(mm2) != "(Intercept)", drop=FALSE]
    levs2 <- levels(ss2[[var]])
    baseline2 <- levs2[1]  # should be g2
    nonbase2 <- levs2[-1]
    level_to_col2 <- setNames(colnames(mm2), nonbase2)
    
    col_drop <- level_to_col2[[g1]]
    if (is.na(col_drop) || col_drop == "") {
      stop(sprintf("Could not find contrast column for %s: %s vs %s.", var, g1, g2))
    }
    if (!col_drop %in% colnames(X_full2)) {
      stop(sprintf("Column to drop '%s' not found in relevelled full matrix.", col_drop))
    }
    
    X_red <- X_full2[, colnames(X_full2) != col_drop, drop=FALSE]
    df_lrt <- ncol(X_full2) - ncol(X_red)
    if (df_lrt != 1) {
      stop(sprintf("Categorical test '%s' did not yield df=1.", test_name))
    }
    
    df_red <- X_red |>
      as.data.frame() |>
      rownames_to_column("sample")
    readr::write_csv(df_red, restricted_path)
    
    group1_res <- g1
    group2_res <- g2
  }
  
  # ---- Store metadata row ----
  lrt_info <- add_row(
    lrt_info,
    test_name         = test_name,
    variable          = var,
    test_type         = test_type,
    group1            = group1_res,
    group2            = group2_res,
    restricted_matrix = file.path("reduced_design_matrices",
                                  paste0("design_matrix_", test_name, ".csv")),
    df_lrt            = 1L,
    lfc_col1          = xcol1,
    lfc_col2          = xcol2
  )
}

readr::write_csv(lrt_info, file.path(outdir, "lrt_tests_metadata.csv"))
