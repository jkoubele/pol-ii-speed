library(rtracklayer)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--gtf",
                    default = '/home/jakub/Desktop/reference_genomes/Ensembl/GRCh38_115/Homo_sapiens.GRCh38.115.gtf',
                     help = "Path to GTF file")
parser$add_argument("--output_dir",
                    default='/home/jakub/Desktop/pol-ii-speed/epigenetics/genomic_intervals',
                    help = "Output directory for BED files")
args <- parser$parse_args()

gtf_path   <- args$gtf
output_dir <- args$output_dir

gtf <- import(gtf_path)

# protein-coding genes only
genes_gr <- gtf[gtf$type == "gene" &
                !is.na(gtf$gene_biotype) &
                gtf$gene_biotype == "protein_coding"]
names(genes_gr) <- genes_gr$gene_id

# strip GTF metadata, keep only what BED export needs
mcols(genes_gr) <- NULL
genes_gr$name  <- names(genes_gr)
genes_gr$score <- 0L

# promoters: symmetric ±2kb around TSS
promoters_gr <- trim(promoters(genes_gr, upstream = 2000, downstream = 2000))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
export(genes_gr,    file.path(output_dir, "gene_bodies.bed"), format = "BED")
export(promoters_gr, file.path(output_dir, "promoters.bed"),  format = "BED")

cat("Done.\n")
cat(sprintf("  gene bodies: %d genes\n", length(genes_gr)))
cat(sprintf("  promoters:   %d regions\n", length(promoters_gr)))