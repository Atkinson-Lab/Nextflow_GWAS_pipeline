#!/usr/bin/env Rscript
#' Create QQ plot from GWAS summary statistics

suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
  library(ggplot2)
  library(argparse)
})

# Parse arguments
parser <- ArgumentParser(description = "Create QQ plot")
parser$add_argument("--sumstats", required = TRUE, help = "Summary statistics file")
parser$add_argument("--lambda-gc", required = FALSE, help = "Lambda GC file")
parser$add_argument("--output-prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--title", default = "QQ Plot", help = "Plot title")
args <- parser$parse_args()

# Read data
cat("Reading summary statistics...\n")
if (grepl("\\.gz$", args$sumstats)) {
  df <- fread(cmd = paste("zcat", args$sumstats))
} else {
  df <- fread(args$sumstats)
}

# Detect P-value column
p_col <- names(df)[grep("^(P|PVALUE|P_VALUE|p.value)$", names(df), ignore.case = TRUE)[1]]

if (is.na(p_col)) {
  stop("Could not detect P-value column")
}

# Extract P-values
p_values <- df[[p_col]]
p_values <- p_values[!is.na(p_values) & p_values > 0 & p_values <= 1]

# Calculate lambda GC
chisq <- qchisq(1 - p_values, 1)
lambda_gc <- median(chisq) / qchisq(0.5, 1)

# Read provided lambda if available
if (!is.null(args$lambda_gc) && file.exists(args$lambda_gc)) {
  lambda_gc <- as.numeric(readLines(args$lambda_gc)[1])
}

# Create title with lambda
plot_title <- paste0(args$title, " (λGC = ", round(lambda_gc, 3), ")")

# Create PNG
cat("Creating QQ plot...\n")
png(paste0(args$output_prefix, ".qq.png"), width = 800, height = 800, res = 100)
qq(p_values, main = plot_title)
dev.off()

# Create PDF
pdf(paste0(args$output_prefix, ".qq.pdf"), width = 8, height = 8)
qq(p_values, main = plot_title)
dev.off()

# Also create a ggplot2 version with confidence intervals
observed <- -log10(sort(p_values))
expected <- -log10(ppoints(length(p_values)))

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# Confidence intervals
n <- length(p_values)
qq_df$lower <- -log10(qbeta(0.975, 1:n, n:1))
qq_df$upper <- -log10(qbeta(0.025, 1:n, n:1))

p <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(size = 0.8, alpha = 0.6, color = "#1565C0") +
  labs(
    title = plot_title,
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(args$output_prefix, ".qq_ggplot.png"), p, width = 8, height = 8, dpi = 100)
ggsave(paste0(args$output_prefix, ".qq_ggplot.pdf"), p, width = 8, height = 8)

cat("QQ plots created successfully\n")
cat("Lambda GC:", round(lambda_gc, 4), "\n")
