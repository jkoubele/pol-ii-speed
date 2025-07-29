#!/usr/bin/env Rscript

library(argparse)
library(rtracklayer)
library(readr)
library(tibble)
library(arrow)


parser <- ArgumentParser()

parser$add_argument("--bed_graph_plus",
                    required = TRUE,
                    help = "A .bedGraph.gz file with coverage for plus strand.")
parser$add_argument("--bed_graph_minus",
                    required = TRUE,
                    help = "A .bedGraph.gz file with coverage for minus strand.")
parser$add_argument("--introns_bed_file",
                    required = TRUE,
                    help = "File with introns in .bed format.")
parser$add_argument("--output_file_basename",
                    required = TRUE,
                    help = "Base name of the output file (without the extension).")
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.",
                    default = '.')

args <- parser$parse_args()


introns <- read_tsv(args$introns_bed_file,
                    col_names = c("chromosome", "start", "end", "name", "score", "strand")
)

bed_graph_plus <- import.bedGraph(args$bed_graph_plus)
bed_graph_minus <- import.bedGraph(args$bed_graph_minus)

coverage_plus <- coverage(bed_graph_plus, weight = bed_graph_plus$score)
coverage_minus <- coverage(bed_graph_minus, weight = bed_graph_minus$score)

output_rescaled_coverages <- list()

num_bins <- 100

for (intron_number in seq_len(nrow(introns))) {
  row <- introns[intron_number,]
  if (intron_number %% 1000 == 0) {
    message(sprintf("[%s] Processed %d introns", Sys.time(), intron_number))
  }

  if (row$strand == "+") {
    coverage_rle <- stats::window(
      coverage_plus[[row$chromosome]],
      start = row$start,
      end = row$end
    )
  } else if (row$strand == "-") {
    coverage_rle <- stats::window(
      coverage_minus[[row$chromosome]],
      start = row$start,
      end = row$end
    )
  } else {
    stop("Invalid strand: ", row$strand)
  }

  coverage <- as.numeric(coverage_rle)
  rescaled_coverage <- approx(
    x = seq_along(coverage),
    y = coverage,
    xout = seq(1, length(coverage), length.out = num_bins)
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
colnames(df_out) <- c("intron_name", sprintf("coverage_bin_%03d", 1:num_bins))
write_parquet(df_out, file.path(args$output_folder, paste0(args$output_file_basename, ".parquet")))
