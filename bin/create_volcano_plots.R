#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("--test_results",
                    required = TRUE,
                    help = "Path to CSV file with LRTs results.")
parser$add_argument("--output_folder", default = '.')

args <- parser$parse_args()

test_results <- read_csv(args$test_results)

min_p_value <- 1e-20
fdr_threshold <- 0.05
l2fc_clamp_value <- 5

for (plot_parameter_type in unique(test_results$parameter_type)) {
  output_subfolder <- file.path(args$output_folder, plot_parameter_type)
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder, recursive = TRUE)
  }
  for (plot_test_id in unique(test_results$test_id)) {
    plot_test_results <- test_results |>
      filter(parameter_type == plot_parameter_type) |>
      filter(test_id == plot_test_id) |>
      mutate(p_value_clipped = pmax(p_value, min_p_value),
             fdr = p.adjust(p_value_clipped, method = 'fdr'),
             category = case_when(
               fdr < fdr_threshold & l2fc_regularized > 0 ~ "up",
               fdr < fdr_threshold & l2fc_regularized < 0 ~ "down",
               TRUE ~ "not_significant"
             ),
             neglog10_fdr = -log10(fdr),
             l2cf_clamped = pmin(pmax(l2fc_regularized, -l2fc_clamp_value), l2fc_clamp_value))

    stopifnot(length(unique(plot_test_results$variable)) == 1)
    stopifnot(length(unique(plot_test_results$group_1)) <= 1)
    stopifnot(length(unique(plot_test_results$group_2)) <= 1)

    cat_counts <- table(plot_test_results$category)
    cat_text <- paste0(
      "down: ", cat_counts["down"], " | ",
      "not significant: ", cat_counts["not_significant"], " | ",
      "up: ", cat_counts["up"]
    )

    group_1 <- plot_test_results |>
      summarize(v = if (all(is.na(group_1))) NA else first(group_1, default = NA)) |>
      pull(v)

    group_2 <- plot_test_results |>
      summarize(v = if (all(is.na(group_2))) NA else first(group_2, default = NA)) |>
      pull(v)

    tested_variable <- plot_test_results |>
      summarize(v = first(variable)) |>
      pull(v)

    title_nice <- str_replace_all(plot_parameter_type, "_", " ") |>
      str_to_sentence() |>
      paste('-', tested_variable)

    if (!is.na(group_1) && !is.na(group_2)) {
      title_nice <- paste0(title_nice, ': ', group_1, ' vs. ', group_2)
    }

    volcano_plot <- ggplot(plot_test_results, aes(x = l2cf_clamped, y = neglog10_fdr, color = category)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(values = c(
        "not_significant" = "gray70",
        "up" = "red",
        "down" = "blue"
      )) +
      labs(
        title = title_nice,
        subtitle = cat_text,
        x = "L2FC estimate",
        y = "-log10(FDR)",
        color = ""
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )

    file_name <- tested_variable
    if (!is.na(group_1) && !is.na(group_2)) {
      file_name <- paste0(plot_parameter_type, '_', file_name, '_', group_1, '_vs_', group_2)
    }

    ggsave(file.path(output_subfolder, paste0(file_name, '.png')),
           plot = volcano_plot,
           bg = 'white',
           width = 10,
           height = 8,
           dpi = 300)

  }
}


