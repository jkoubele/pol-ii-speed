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
parser$add_argument("--introns_count_file",
                    required = TRUE,
                    help = "File with intron reads counts.")
parser$add_argument("--output_file_basename",
                    required = TRUE,
                    help = "Base name of the output file (without the extension).")
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.",
                    default = '.')
parser$add_argument("--num_bins",
                    default = 100,
                    type = "integer",
                    help = "Number of bins into which the coverage will be rescaled. ")

args <- parser$parse_args()

introns <- read_tsv(args$introns_bed_file,
                    col_names = c("chromosome", "start", "end", "name", "score", "strand")
)

bed_graph_plus <- import.bedGraph(args$bed_graph_plus)
bed_graph_minus <- import.bedGraph(args$bed_graph_minus)

num_bins <- args$num_bins

coverage_plus <- coverage(bed_graph_plus, weight = bed_graph_plus$score)
coverage_minus <- coverage(bed_graph_minus, weight = bed_graph_minus$score)

intron_count <- read_tsv(
  args$introns_count_file,
  show_col_types = FALSE
)
intron_count_lookup <- setNames(intron_count$count, intron_count$name)

output_rescaled_coverages <- matrix(
  0,
  nrow = nrow(introns),
  ncol = num_bins
)
rownames(output_rescaled_coverages) <- introns$name

for (intron_number in seq_len(nrow(introns))) {

  chromosome <- introns$chromosome[intron_number]
  start <- introns$start[intron_number]
  end <- introns$end[intron_number]
  strand <- introns$strand[intron_number]
  intron_name <- introns$name[intron_number]

  if (intron_number %% 1000 == 0) {
    message(sprintf("[%s] Processed %d introns", Sys.time(), intron_number))
  }

  if (strand == "+") {
    coverage_rle <- stats::window(
      coverage_plus[[chromosome]],
      start = start + 1, # Adding +1 to start since BED is 0-indexed (half-open) and window() is 1-indexed (closed interval)
      end = end
    )
  } else if (strand == "-") {
    coverage_rle <- stats::window(
      coverage_minus[[chromosome]],
      start = start + 1, # Adding +1 to start since BED is 0-indexed (half-open) and window() is 1-indexed (closed interval)
      end = end
    )
  } else {
    stop("Invalid strand: ", strand)
  }


  coverage <- as.numeric(coverage_rle)
  rescaled_coverage <- approx(
    x = seq_along(coverage),
    y = coverage,
    xout = seq(1, length(coverage), length.out = num_bins)
  )$y

  rescaled_sum <- sum(rescaled_coverage)

  if (rescaled_sum > 0) {
    rescaled_coverage <- rescaled_coverage * intron_count_lookup[[intron_name]] / rescaled_sum
  }

  if (strand == "-") {
    rescaled_coverage <- rev(rescaled_coverage)
  }

  output_rescaled_coverages[intron_number,] <- rescaled_coverage

}

colnames(output_rescaled_coverages) <- sprintf("coverage_bin_%03d", 1:num_bins)
df_out <- as_tibble(output_rescaled_coverages, rownames = "intron_name")
write_parquet(df_out, file.path(args$output_folder, paste0(args$output_file_basename, ".parquet")))
