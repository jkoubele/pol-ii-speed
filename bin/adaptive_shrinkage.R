#!/usr/bin/env Rscript

library(tidyverse)
library(ashr)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("--model_parameters",
                    required = TRUE,
                    help = "Path to CSV file with model parameters.")
parser$add_argument("--output_folder", default = ".",
                    help = "Where to write the output CSV with regularization coefficients.")

# if (interactive()) {
#   args = list(model_parameters = '/cellfile/projects/pol_ii_speed/jkoubele/analysis/mouse_age_dr/modeling/model_run_02_Dec_2025_18_55_48/model_results/model_parameters.csv',
#               output_folder = '/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/test_ashr')
# } else {
#   args <- parser$parse_args()
# }

args <- parser$parse_args()

model_parameters <- read.csv(args$model_parameters)

regularization_coefficients_df <- tibble(
  parameter_type = character(),
  feature_name = character(),
  prior_sd = double(),
  lambda = double()
)

for (parameter in c('beta', 'gamma')) {
  df_parameter <- model_parameters |> filter(parameter_type == parameter)
  for (feature in unique(df_parameter$feature_name)) {
    df_feature <- df_parameter |>
      filter(feature_name == feature,
             identifiable == TRUE,
             training_diverged_full_model == FALSE,
             !is.na(value),
             !is.na(SE),
             SE > 0)

    if (nrow(df_feature) == 0) {
      stop(sprintf(
        "No parameter LFC and SE available for parameter_type='%s', feature_name='%s'; cannot perform adaptive shrinkage.",
        parameter,
        feature
      ))
    }

    ash_results <- ash(betahat = df_feature$value, sebetahat = df_feature$SE)
    prior_sd <- calc_mixsd(ash_results$fitted_g)

    regularization_coefficients_df <- regularization_coefficients_df |> add_row(
      parameter_type = parameter,
      feature_name = feature,
      prior_sd = prior_sd,
      lambda = 1 / (2 * prior_sd^2)
    )

  }
}

write_csv(regularization_coefficients_df, file.path(args$output_folder, "regularization_coefficients.csv"))

