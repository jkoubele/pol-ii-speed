library(argparse)
library(rtracklayer)
library(tidyverse)
library(arrow)


parser <- ArgumentParser()
parser$add_argument("--bed_graph_plus",
                    help = "bedGraph file with coverage for plus strand.",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intron_coverage/test/coverage_plus_strand.bedGraph.gz')
parser$add_argument("--bed_graph_minus",
                    help = "bedGraph file with coverage for minus strand.",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intron_coverage/test/coverage_minus_strand.bedGraph.gz')
parser$add_argument("--introns_file",
                    help = "File with introns in .bed format.",
                    default = '/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/introns_filtered.bed')
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.",
                    default = '/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intron_coverage/test')
args <- parser$parse_args()


introns <- read_tsv(args$introns_file,
                    col_names = c("chromosome", "start", "end", "name", "score", "strand")
)
bed_graph_plus <- import.bedGraph(args$bed_graph_plus)
bed_graph_minus <- import.bedGraph(args$bed_graph_minus)
padding <- 10

coverage_plus <- coverage(bed_graph_plus, weight = bed_graph_plus$score)
coverage_minus <- coverage(bed_graph_minus, weight = bed_graph_minus$score)

output_rescaled_coverages <- list()


for (i in seq_len(nrow(introns))) {
  row <- introns[i,]
  if (i %% 1000 == 0) {
    message("Processing coverage of introns: ", i)
  }

  if ((row$end - row$start) <= 2 * padding) {
    next
  }

  if (row$strand == "+") {
    coverage_rle <- stats::window(
      coverage_plus[[row$chromosome]],
      start = row$start + padding,
      end = row$end - padding
    )
  } else if (row$strand == "-") {
    coverage_rle <- stats::window(
      coverage_minus[[row$chromosome]],
      start = row$start + padding,
      end = row$end - padding
    )
  } else {
    stop("Invalid strand: ", row$strand)
  }

  coverage <- as.numeric(coverage_rle)
  rescaled_coverage <- approx(
    x = seq_along(coverage),
    y = coverage,
    xout = seq(1, length(coverage), length.out = 100)
  )$y

  if (sum(rescaled_coverage) > 0) {
    rescaled_coverage <- rescaled_coverage * sum(coverage) / sum(rescaled_coverage)
  }

  if (row$strand == "-") {
    rescaled_coverage <- rev(rescaled_coverage)
  }

  output_rescaled_coverages[[row$name]] <- rescaled_coverage
}

df_out <- as_tibble(do.call(rbind, output_rescaled_coverages), rownames = "intron_name")
write_parquet(df_out, file.path(args$output_folder, "rescaled_intron_coverages.parquet"))

