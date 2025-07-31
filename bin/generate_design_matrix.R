library(argparse)
library(readr)

parser <- ArgumentParser()
parser$add_argument("--samplesheet",
                    help = "A .csv file with sample annotation.",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/sample_annotation/sample_annotation.tsv')
parser$add_argument("--formula",
                    required=TRUE,
                    help = "Design formula using standard R syntax (e.g. '~ age + genotype' or '~ age:condition')")
parser$add_argument("--output_folder",
                    help = "",
                    default = '.')

args <- parser$parse_args()

data <- readr::read_csv(args$samplesheet)

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


