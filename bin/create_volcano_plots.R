library(tidyverse)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("--model_results",
                    required = TRUE,
                    help = "Path to model_results.csv")
parser$add_argument("--output_folder", default = '.')

args <- parser$parse_args()

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

df_results <- read_csv(args$model_results) |> 
  filter(!if_any(c(log2_fold_change, p_value), ~ is.na(.) | is.nan(.) | is.infinite(.)))

min_p_value <- 1e-25
fdr_threshold <- 0.05
lfc_clamp_value <- 5

for (param in c("elongation_speed", "splicing_speed")) {
  for (feature in unique(df_results$feature_name)){
    df_parameter <- df_results |>
      filter(parameter_type == param) |>
      filter(feature_name == feature) |>
      mutate(p_value_clipped =  pmax(p_value, min_p_value),
             fdr = p.adjust(p_value_clipped, method='fdr'),
             category = case_when(
               fdr < fdr_threshold & log2_fold_change > 0 ~ "up",
               fdr < fdr_threshold & log2_fold_change < 0 ~ "down",
               TRUE ~ "not_significant"
             ),
             neglog10_fdr = -log10(fdr),
             lcf_clamped = pmin(pmax(log2_fold_change, -lfc_clamp_value), lfc_clamp_value))
    
    cat_counts <- table(df_parameter$category)
    cat_text <- paste0(
      "down: ", cat_counts["down"], " | ",
      "not significant: ", cat_counts["not_significant"], " | ",
      "up: ", cat_counts["up"]
    )
    
    title_nice <- str_replace_all(param, "_", " ") |>
      str_to_sentence()
    
    volcano_plot <- ggplot(df_parameter, aes(x = lcf_clamped, y = neglog10_fdr, color = category)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(values = c(
        "not_significant" = "gray70",
        "up" = "red",
        "down" = "blue"
      )) +
      labs(
        title = paste(title_nice, " - ", feature),
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
    
    ggsave(file.path(output_folder, paste0(param, '_', feature, '.png')),
           plot = volcano_plot,
           bg = 'white',
           width = 10, 
           height = 8, 
           dpi = 300)
  }
   
}

