#!/usr/bin/env Rscript
#' Validate PRS performance in target cohort

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
  library(argparse)
})

# Parse arguments
parser <- ArgumentParser(description = "Validate PRS performance")
parser$add_argument("--prs-scores", required = TRUE, help = "PRS scores file")
parser$add_argument("--phenotype", required = TRUE, help = "Phenotype file")
parser$add_argument("--trait", required = TRUE, help = "Trait name")
parser$add_argument("--output-prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--binary", default = "false", help = "Binary trait (true/false)")
args <- parser$parse_args()

is_binary <- tolower(args$binary) == "true"

cat("Loading data...\n")

# Read PRS scores
prs <- fread(args$prs_scores)
pheno <- fread(args$phenotype)

# Detect ID column
id_cols <- c("IID", "FID", "sample_id", "SAMPLE_ID", "ID")
prs_id <- names(prs)[names(prs) %in% id_cols][1]
pheno_id <- names(pheno)[names(pheno) %in% id_cols][1]

# Detect PRS column
prs_col <- names(prs)[grep("PRS|SCORE|score", names(prs), ignore.case = TRUE)][1]
if (is.na(prs_col)) {
  prs_col <- names(prs)[ncol(prs)]  # Assume last column is PRS
}

# Merge data
merged <- merge(prs, pheno, by.x = prs_id, by.y = pheno_id)

if (nrow(merged) == 0) {
  stop("No samples matched between PRS and phenotype files")
}

cat("Matched", nrow(merged), "samples\n")

# Get phenotype column
if (!args$trait %in% names(merged)) {
  stop("Trait ", args$trait, " not found in phenotype file")
}

# Prepare analysis data
analysis_df <- data.frame(
  id = merged[[prs_id]],
  prs = merged[[prs_col]],
  pheno = merged[[args$trait]]
)
analysis_df <- analysis_df[complete.cases(analysis_df), ]

cat("Analyzing", nrow(analysis_df), "samples with complete data\n")

# Calculate performance metrics
results <- list()

if (is_binary) {
  # Binary trait: calculate AUC
  roc_obj <- roc(analysis_df$pheno, analysis_df$prs, quiet = TRUE)
  results$auc <- auc(roc_obj)
  results$auc_ci <- ci.auc(roc_obj)

  # Calculate Nagelkerke R²
  null_model <- glm(pheno ~ 1, data = analysis_df, family = binomial)
  full_model <- glm(pheno ~ prs, data = analysis_df, family = binomial)

  n <- nrow(analysis_df)
  results$r2_nagelkerke <- (1 - exp((logLik(null_model) - logLik(full_model))[1] * 2/n)) /
                           (1 - exp(logLik(null_model)[1] * 2/n))

  # Odds ratio per SD
  analysis_df$prs_scaled <- scale(analysis_df$prs)
  or_model <- glm(pheno ~ prs_scaled, data = analysis_df, family = binomial)
  results$or_per_sd <- exp(coef(or_model)[2])
  results$or_ci <- exp(confint(or_model)[2, ])

} else {
  # Quantitative trait: calculate R²
  lm_model <- lm(pheno ~ prs, data = analysis_df)
  results$r2 <- summary(lm_model)$r.squared
  results$adj_r2 <- summary(lm_model)$adj.r.squared

  # Beta coefficient
  results$beta <- coef(lm_model)[2]
  results$beta_se <- summary(lm_model)$coefficients[2, 2]
  results$beta_p <- summary(lm_model)$coefficients[2, 4]

  # Correlation
  results$correlation <- cor(analysis_df$prs, analysis_df$pheno)
}

# Write results
cat("\n=== Results ===\n")

results_df <- data.frame(
  metric = names(results),
  value = sapply(results, function(x) {
    if (length(x) > 1) paste(round(x, 4), collapse = ", ")
    else round(x, 4)
  })
)

print(results_df)

fwrite(results_df, paste0(args$output_prefix, ".validation_results.tsv"), sep = "\t")

# Save as JSON
json_out <- toJSON(results, auto_unbox = TRUE, pretty = TRUE)
writeLines(json_out, paste0(args$output_prefix, ".validation_metrics.json"))

# Create visualization
cat("\nCreating plots...\n")

if (is_binary) {
  # ROC curve
  pdf(paste0(args$output_prefix, ".validation_plot.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))

  # ROC
  plot(roc_obj, main = paste0("ROC Curve (AUC = ", round(results$auc, 3), ")"))

  # PRS distribution by case/control
  boxplot(prs ~ pheno, data = analysis_df,
          main = "PRS by Case/Control Status",
          xlab = "Status", ylab = "PRS",
          col = c("#90CAF9", "#EF9A9A"))

  dev.off()

} else {
  # Scatter plot
  p <- ggplot(analysis_df, aes(x = prs, y = pheno)) +
    geom_point(alpha = 0.5, color = "#1565C0") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(
      title = paste0("PRS vs ", args$trait, " (R² = ", round(results$r2, 3), ")"),
      x = "Polygenic Risk Score",
      y = args$trait
    ) +
    theme_bw()

  ggsave(paste0(args$output_prefix, ".validation_plot.pdf"), p, width = 8, height = 6)
}

cat("Validation complete\n")
