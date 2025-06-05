library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("--sample_annotation",
                    help = ".tsv file with sample annotation",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/sample_annotation/sample_annotation.tsv')
parser$add_argument("--output_folder",
                    help = "",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/sample_annotation')
args <- parser$parse_args()

data <- read_tsv(args$sample_annotation)

# data <- data |> mutate(genotype = fct_relevel(as.factor(genotype), "wt"),
#                        age = fct_relevel(as.factor(age), "young"))
data <- data |> mutate(age = fct_relevel(as.factor(age_group), "young"))

formula <- ~age
design_matrix <- model.matrix(formula, data)
design_matrix <- design_matrix[, colnames(design_matrix) != "(Intercept)", drop = FALSE]

rownames(design_matrix) <- data$sample_name

df_out <- design_matrix |>
  as.data.frame() |>
  rownames_to_column(var = "sample_name")

write_tsv(df_out, file.path(args$output_folder, 'design_matrix.tsv'))


