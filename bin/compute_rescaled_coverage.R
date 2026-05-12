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
parser$add_argument("--edge_margin",
                    default = 10,
                    type = "integer",
                    help = paste("Number of bp trimmed from each intron edge before resampling.",
                                 "Should match the --overlap_bp_threshold used in extract_intronic_reads.py",
                                 "so the trimmed region contains only positions where intron-classified reads",
                                 "are sampled without boundary bias."))

args <- parser$parse_args()

introns <- read_tsv(args$introns_bed_file,
                    col_names = c("chromosome", "start", "end", "name", "score", "strand"),
                    show_col_types = FALSE
)

bed_graph_plus <- import.bedGraph(args$bed_graph_plus)
bed_graph_minus <- import.bedGraph(args$bed_graph_minus)

num_bins <- args$num_bins
edge_margin <- args$edge_margin

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

  # Trim `edge_margin` bp from each side before resampling. BED is 0-indexed
  # half-open, window() is 1-indexed closed, so the full intron is
  # window(start + 1, end) and the trimmed region is window(start + 1 + edge_margin, end - edge_margin).
  trim_start <- start + 1 + edge_margin
  trim_end <- end - edge_margin

  if (trim_end - trim_start < 1) {
    next
  }

  if (strand == "+") {
    coverage_rle <- stats::window(
      coverage_plus[[chromosome]],
      start = trim_start,
      end = trim_end
    )
  } else if (strand == "-") {
    coverage_rle <- stats::window(
      coverage_minus[[chromosome]],
      start = trim_start,
      end = trim_end
    )
  } else {
    stop("Invalid strand: ", strand)
  }


  coverage <- as.numeric(coverage_rle)
  L <- length(coverage)
  if (L <= num_bins) {
    rescaled_coverage <- approx(
      x = seq_along(coverage),
      y = coverage,
      xout = seq(1, L, length.out = num_bins)
    )$y
  } else {
    bin_edges <- round(seq(0L, L, length.out = num_bins + 1L))
    rescaled_coverage <- vapply(seq_len(num_bins), function(b) {
      mean(coverage[(bin_edges[b] + 1L):bin_edges[b + 1L]])
    }, FUN.VALUE = numeric(1))
  }

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
