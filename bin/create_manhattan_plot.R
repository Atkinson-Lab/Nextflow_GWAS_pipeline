#!/usr/bin/env Rscript
#' Create Manhattan plot from GWAS summary statistics

suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
  library(ggplot2)
  library(argparse)
})

# Parse arguments
parser <- ArgumentParser(description = "Create Manhattan plot")
parser$add_argument("--sumstats", required = TRUE, help = "Summary statistics file")
parser$add_argument("--lambda-gc", required = FALSE, help = "Lambda GC file")
parser$add_argument("--output-prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--title", default = "Manhattan Plot", help = "Plot title")
parser$add_argument("--genome-wide", type = "double", default = 5e-8, help = "Genome-wide significance threshold")
parser$add_argument("--suggestive", type = "double", default = 1e-5, help = "Suggestive threshold")
args <- parser$parse_args()

# Read data
cat("Reading summary statistics...\n")
if (grepl("\\.gz$", args$sumstats)) {
  df <- fread(cmd = paste("zcat", args$sumstats))
} else {
  df <- fread(args$sumstats)
}

# Detect columns
chr_col <- names(df)[grep("^(CHR|CHROM|#CHROM|chromosome)$", names(df), ignore.case = TRUE)[1]]
pos_col <- names(df)[grep("^(POS|BP|GENPOS|position)$", names(df), ignore.case = TRUE)[1]]
snp_col <- names(df)[grep("^(SNP|ID|RSID|variant_id|MarkerName)$", names(df), ignore.case = TRUE)[1]]
p_col <- names(df)[grep("^(P|PVALUE|P_VALUE|p.value)$", names(df), ignore.case = TRUE)[1]]

if (is.na(chr_col) || is.na(pos_col) || is.na(p_col)) {
  stop("Could not detect required columns (CHR, POS, P)")
}

# Prepare data for manhattan plot
plot_df <- data.frame(
  CHR = df[[chr_col]],
  BP = df[[pos_col]],
  P = df[[p_col]],
  SNP = if (!is.na(snp_col)) df[[snp_col]] else paste0("chr", df[[chr_col]], ":", df[[pos_col]])
)

# Remove invalid entries
plot_df <- plot_df[!is.na(plot_df$P) & plot_df$P > 0 & plot_df$P <= 1, ]
plot_df <- plot_df[!is.na(plot_df$CHR) & !is.na(plot_df$BP), ]

# Convert chromosome to numeric
plot_df$CHR <- as.numeric(gsub("chr", "", plot_df$CHR, ignore.case = TRUE))
plot_df <- plot_df[!is.na(plot_df$CHR), ]

# Read lambda GC if provided
lambda_gc <- NA
if (!is.null(args$lambda_gc) && file.exists(args$lambda_gc)) {
  lambda_gc <- as.numeric(readLines(args$lambda_gc)[1])
}

# Create title with lambda
plot_title <- args$title
if (!is.na(lambda_gc)) {
  plot_title <- paste0(plot_title, " (λGC = ", round(lambda_gc, 3), ")")
}

# Create PNG
cat("Creating Manhattan plot...\n")
png(paste0(args$output_prefix, ".manhattan.png"), width = 1200, height = 600, res = 100)
manhattan(
  plot_df,
  main = plot_title,
  genomewideline = -log10(args$genome_wide),
  suggestiveline = -log10(args$suggestive),
  col = c("#1565C0", "#90CAF9"),
  cex = 0.6,
  cex.axis = 0.9
)
dev.off()

# Create PDF
pdf(paste0(args$output_prefix, ".manhattan.pdf"), width = 12, height = 6)
manhattan(
  plot_df,
  main = plot_title,
  genomewideline = -log10(args$genome_wide),
  suggestiveline = -log10(args$suggestive),
  col = c("#1565C0", "#90CAF9"),
  cex = 0.6,
  cex.axis = 0.9
)
dev.off()

cat("Manhattan plot created successfully\n")
