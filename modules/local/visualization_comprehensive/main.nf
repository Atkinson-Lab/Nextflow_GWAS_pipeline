/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPREHENSIVE VISUALIZATION MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Publication-quality visualizations for:
    - Colocalization analysis (heatmaps, locus plots, PP4 summaries)
    - Cross-ancestry heritability comparison
    - PRS accuracy by ancestry and tool
    - Multi-ancestry GWAS comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COLOCALIZATION_PLOTS {
    tag "$meta.trait"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2=3.4 conda-forge::r-pheatmap conda-forge::r-viridis conda-forge::r-cowplot conda-forge::r-dplyr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    tuple val(meta), path(coloc_results)
    val qtl_types     // List of QTL types analyzed

    output:
    tuple val(meta), path("${prefix}.coloc_heatmap.pdf"), emit: heatmap
    tuple val(meta), path("${prefix}.coloc_pp4_barplot.pdf"), emit: pp4_barplot
    tuple val(meta), path("${prefix}.coloc_tissue_summary.pdf"), emit: tissue_summary
    tuple val(meta), path("${prefix}.coloc_locuszoom.pdf"), emit: locuszoom, optional: true
    tuple val(meta), path("${prefix}.coloc_summary.tsv"), emit: summary_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)
    library(pheatmap)
    library(viridis)
    library(cowplot)
    library(tidyr)

    # Read colocalization results
    coloc <- read.table("${coloc_results}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

    # =========================================================================
    # 1. PP4 HEATMAP: Genes vs QTL Types/Tissues
    # =========================================================================
    if (nrow(coloc) > 0 && "gene" %in% names(coloc) && "qtl_type" %in% names(coloc)) {
        # Create matrix for heatmap
        coloc_wide <- coloc %>%
            select(gene, qtl_type, PP.H4.abf) %>%
            pivot_wider(names_from = qtl_type, values_from = PP.H4.abf, values_fn = max)

        mat <- as.matrix(coloc_wide[,-1])
        rownames(mat) <- coloc_wide\$gene

        # Replace NA with 0
        mat[is.na(mat)] <- 0

        # Only plot if we have data
        if (nrow(mat) > 1 && ncol(mat) > 1) {
            pdf("${prefix}.coloc_heatmap.pdf", width=12, height=max(8, nrow(mat)*0.3))
            pheatmap(mat,
                     main = "Colocalization PP4 by Gene and QTL Type\\n${meta.trait}",
                     color = viridis(100),
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     fontsize_row = 8,
                     fontsize_col = 10,
                     display_numbers = TRUE,
                     number_format = "%.2f",
                     number_color = "white",
                     annotation_legend = TRUE)
            dev.off()
        } else {
            # Placeholder
            pdf("${prefix}.coloc_heatmap.pdf", width=8, height=6)
            plot.new()
            text(0.5, 0.5, "Insufficient data for heatmap", cex=1.5)
            dev.off()
        }
    }

    # =========================================================================
    # 2. PP4 BARPLOT: Top colocalizing genes
    # =========================================================================
    top_coloc <- coloc %>%
        filter(PP.H4.abf >= 0.5) %>%
        arrange(desc(PP.H4.abf)) %>%
        head(30)

    if (nrow(top_coloc) > 0) {
        p_bar <- ggplot(top_coloc, aes(x = reorder(gene, PP.H4.abf), y = PP.H4.abf, fill = qtl_type)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_viridis_d() +
            labs(title = "Top Colocalizing Genes (PP4 >= 0.5)",
                 subtitle = "${meta.trait}",
                 x = "Gene",
                 y = "Posterior Probability of Colocalization (PP4)",
                 fill = "QTL Type") +
            theme_minimal() +
            theme(legend.position = "right",
                  axis.text.y = element_text(size = 8))

        ggsave("${prefix}.coloc_pp4_barplot.pdf", p_bar, width = 10, height = max(6, nrow(top_coloc)*0.25))
    } else {
        pdf("${prefix}.coloc_pp4_barplot.pdf", width=8, height=6)
        plot.new()
        text(0.5, 0.5, "No colocalizations with PP4 >= 0.5", cex=1.5)
        dev.off()
    }

    # =========================================================================
    # 3. TISSUE/CELL TYPE SUMMARY
    # =========================================================================
    if ("tissue" %in% names(coloc)) {
        tissue_summary <- coloc %>%
            filter(PP.H4.abf >= 0.5) %>%
            group_by(tissue, qtl_type) %>%
            summarize(n_coloc = n(), mean_pp4 = mean(PP.H4.abf), .groups = "drop")

        if (nrow(tissue_summary) > 0) {
            p_tissue <- ggplot(tissue_summary, aes(x = tissue, y = n_coloc, fill = qtl_type)) +
                geom_bar(stat = "identity", position = "dodge") +
                coord_flip() +
                scale_fill_viridis_d() +
                labs(title = "Colocalizations by Tissue/Cell Type",
                     subtitle = "${meta.trait} (PP4 >= 0.5)",
                     x = "Tissue/Cell Type",
                     y = "Number of Colocalizations",
                     fill = "QTL Type") +
                theme_minimal()

            ggsave("${prefix}.coloc_tissue_summary.pdf", p_tissue, width = 10, height = 8)
        }
    } else {
        pdf("${prefix}.coloc_tissue_summary.pdf", width=8, height=6)
        plot.new()
        text(0.5, 0.5, "No tissue information available", cex=1.5)
        dev.off()
    }

    # =========================================================================
    # 4. SUMMARY TABLE
    # =========================================================================
    summary_stats <- coloc %>%
        summarize(
            n_loci_tested = n(),
            n_coloc_pp4_50 = sum(PP.H4.abf >= 0.5),
            n_coloc_pp4_75 = sum(PP.H4.abf >= 0.75),
            n_coloc_pp4_90 = sum(PP.H4.abf >= 0.90),
            mean_pp4 = mean(PP.H4.abf),
            max_pp4 = max(PP.H4.abf)
        )

    write.table(summary_stats, "${prefix}.coloc_summary.tsv", sep="\\t", row.names=FALSE, quote=FALSE)

    # Version info
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    ggplot2: "', packageVersion("ggplot2"), '"')
    ), "versions.yml")
    """
}

process HERITABILITY_COMPARISON_PLOTS {
    tag "h2_comparison"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2=3.4 conda-forge::r-dplyr conda-forge::r-tidyr conda-forge::r-scales"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    path h2_results       // All heritability results combined
    val traits            // List of traits

    output:
    path "h2_cross_ancestry_comparison.pdf", emit: h2_comparison
    path "h2_forest_plot.pdf", emit: forest_plot
    path "h2_enrichment_plot.pdf", emit: enrichment, optional: true
    path "rg_cross_ancestry_heatmap.pdf", emit: rg_heatmap, optional: true
    path "h2_summary_table.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)

    # Read heritability results
    h2 <- read.table("${h2_results}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

    # Define ancestry colors
    ancestry_colors <- c(
        "EUR" = "#4DAF4A",
        "AFR" = "#FF7F00",
        "EAS" = "#377EB8",
        "SAS" = "#984EA3",
        "AMR" = "#E41A1C",
        "MID" = "#A65628",
        "AAC" = "#F781BF",
        "AHI" = "#999999",
        "LATINO" = "#FFFF33",
        "META" = "#000000"
    )

    # =========================================================================
    # 1. CROSS-ANCESTRY H2 COMPARISON (Bar Plot with Error Bars)
    # =========================================================================
    p_h2 <- ggplot(h2, aes(x = ancestry, y = h2, fill = ancestry)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        geom_errorbar(aes(ymin = h2 - h2_se, ymax = h2 + h2_se),
                      width = 0.2, position = position_dodge(0.7)) +
        facet_wrap(~trait, scales = "free_y", ncol = 3) +
        scale_fill_manual(values = ancestry_colors) +
        labs(title = "SNP-Heritability by Ancestry Group",
             subtitle = "LDSC estimates with standard errors",
             x = "Ancestry",
             y = expression(paste("SNP-heritability (", h^2, ")")),
             fill = "Ancestry") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom",
              strip.background = element_rect(fill = "lightblue"))

    ggsave("h2_cross_ancestry_comparison.pdf", p_h2, width = 12, height = 8)

    # =========================================================================
    # 2. FOREST PLOT
    # =========================================================================
    h2_ordered <- h2 %>%
        arrange(trait, h2) %>%
        mutate(label = paste(trait, ancestry, sep = " - "))

    p_forest <- ggplot(h2_ordered, aes(x = h2, y = reorder(label, h2), color = ancestry)) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = h2 - 1.96*h2_se, xmax = h2 + 1.96*h2_se), height = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        scale_color_manual(values = ancestry_colors) +
        labs(title = "Heritability Forest Plot",
             subtitle = "Point estimates with 95% confidence intervals",
             x = expression(paste("SNP-heritability (", h^2, ")")),
             y = "",
             color = "Ancestry") +
        theme_bw() +
        theme(legend.position = "right",
              axis.text.y = element_text(size = 8))

    ggsave("h2_forest_plot.pdf", p_forest, width = 10, height = max(6, nrow(h2)*0.3))

    # =========================================================================
    # 3. GENETIC CORRELATION HEATMAP (if rg data available)
    # =========================================================================
    if (file.exists("rg_results.tsv")) {
        rg <- read.table("rg_results.tsv", header=TRUE, sep="\\t")

        rg_wide <- rg %>%
            select(ancestry1, ancestry2, rg, trait) %>%
            pivot_wider(names_from = ancestry2, values_from = rg)

        # Create separate heatmaps per trait
        pdf("rg_cross_ancestry_heatmap.pdf", width = 10, height = 8)
        for (t in unique(rg\$trait)) {
            rg_t <- rg %>% filter(trait == t)
            mat <- matrix(NA, length(unique(rg_t\$ancestry1)), length(unique(rg_t\$ancestry2)))
            rownames(mat) <- unique(rg_t\$ancestry1)
            colnames(mat) <- unique(rg_t\$ancestry2)
            for (i in 1:nrow(rg_t)) {
                mat[rg_t\$ancestry1[i], rg_t\$ancestry2[i]] <- rg_t\$rg[i]
            }
            heatmap(mat, main = paste("Cross-Ancestry Genetic Correlation:", t),
                    col = colorRampPalette(c("blue", "white", "red"))(100),
                    scale = "none")
        }
        dev.off()
    }

    # =========================================================================
    # 4. SUMMARY TABLE
    # =========================================================================
    summary_h2 <- h2 %>%
        group_by(trait) %>%
        summarize(
            n_ancestries = n(),
            mean_h2 = mean(h2, na.rm=TRUE),
            sd_h2 = sd(h2, na.rm=TRUE),
            max_h2 = max(h2, na.rm=TRUE),
            max_h2_ancestry = ancestry[which.max(h2)],
            min_h2 = min(h2, na.rm=TRUE),
            min_h2_ancestry = ancestry[which.min(h2)],
            h2_range = max_h2 - min_h2
        )

    write.table(summary_h2, "h2_summary_table.tsv", sep="\\t", row.names=FALSE, quote=FALSE)

    # Version info
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    ggplot2: "', packageVersion("ggplot2"), '"')
    ), "versions.yml")
    """
}

process PRS_COMPARISON_PLOTS {
    tag "prs_comparison"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2=3.4 conda-forge::r-dplyr conda-forge::r-tidyr conda-forge::r-patchwork"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    path prs_validation_results    // Combined PRS validation results

    output:
    path "prs_r2_by_ancestry.pdf", emit: r2_ancestry
    path "prs_r2_by_method.pdf", emit: r2_method
    path "prs_method_ancestry_heatmap.pdf", emit: heatmap
    path "prs_auc_comparison.pdf", emit: auc, optional: true
    path "prs_calibration.pdf", emit: calibration, optional: true
    path "prs_best_method_summary.pdf", emit: best_method
    path "prs_correlation_matrix.pdf", emit: correlation
    path "prs_summary_table.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)

    # Read PRS validation results
    prs <- read.table("${prs_validation_results}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

    # Define colors
    method_colors <- c(
        "PRS-CSx" = "#E41A1C",
        "PRS-CS" = "#377EB8",
        "GAUDI" = "#4DAF4A",
        "LDpred2" = "#984EA3",
        "PRSice2" = "#FF7F00",
        "C+T" = "#FFFF33"
    )

    ancestry_colors <- c(
        "EUR" = "#4DAF4A",
        "AFR" = "#FF7F00",
        "EAS" = "#377EB8",
        "SAS" = "#984EA3",
        "AMR" = "#E41A1C",
        "AAC" = "#F781BF",
        "LATINO" = "#999999"
    )

    # =========================================================================
    # 1. R² BY ANCESTRY (Grouped Bar Plot)
    # =========================================================================
    p_r2_ancestry <- ggplot(prs, aes(x = ancestry, y = r2, fill = method)) +
        geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
        geom_errorbar(aes(ymin = r2 - r2_se, ymax = r2 + r2_se),
                      position = position_dodge(0.8), width = 0.2) +
        facet_wrap(~trait, scales = "free_y", ncol = 2) +
        scale_fill_manual(values = method_colors) +
        labs(title = "PRS Prediction Accuracy by Ancestry",
             subtitle = "Incremental R² with standard errors",
             x = "Ancestry Group",
             y = expression(paste("Incremental ", R^2)),
             fill = "PRS Method") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom",
              strip.background = element_rect(fill = "lightblue"))

    ggsave("prs_r2_by_ancestry.pdf", p_r2_ancestry, width = 14, height = 10)

    # =========================================================================
    # 2. R² BY METHOD (Shows performance gap across ancestries)
    # =========================================================================
    p_r2_method <- ggplot(prs, aes(x = method, y = r2, fill = ancestry)) +
        geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
        geom_errorbar(aes(ymin = r2 - r2_se, ymax = r2 + r2_se),
                      position = position_dodge(0.8), width = 0.2) +
        facet_wrap(~trait, scales = "free_y", ncol = 2) +
        scale_fill_manual(values = ancestry_colors) +
        labs(title = "PRS Method Performance Across Ancestries",
             subtitle = "Incremental R² - Shows transferability gap",
             x = "PRS Method",
             y = expression(paste("Incremental ", R^2)),
             fill = "Ancestry") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")

    ggsave("prs_r2_by_method.pdf", p_r2_method, width = 14, height = 10)

    # =========================================================================
    # 3. METHOD x ANCESTRY HEATMAP
    # =========================================================================
    prs_wide <- prs %>%
        select(method, ancestry, r2, trait) %>%
        pivot_wider(names_from = ancestry, values_from = r2)

    # Create heatmap for each trait
    pdf("prs_method_ancestry_heatmap.pdf", width = 12, height = 8)
    for (t in unique(prs\$trait)) {
        prs_t <- prs %>% filter(trait == t)

        p_heat <- ggplot(prs_t, aes(x = ancestry, y = method, fill = r2)) +
            geom_tile(color = "white") +
            geom_text(aes(label = sprintf("%.3f", r2)), color = "black", size = 3) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                 midpoint = median(prs_t\$r2, na.rm=TRUE)) +
            labs(title = paste("PRS Prediction Accuracy:", t),
                 subtitle = "Incremental R² by Method and Ancestry",
                 x = "Ancestry", y = "PRS Method",
                 fill = expression(R^2)) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        print(p_heat)
    }
    dev.off()

    # =========================================================================
    # 4. BEST METHOD BY ANCESTRY
    # =========================================================================
    best_method <- prs %>%
        group_by(trait, ancestry) %>%
        slice_max(r2, n = 1) %>%
        ungroup()

    p_best <- ggplot(best_method, aes(x = ancestry, y = trait, fill = method)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.3f", r2)), color = "black", size = 3) +
        scale_fill_manual(values = method_colors) +
        labs(title = "Best PRS Method by Ancestry and Trait",
             subtitle = "Method with highest R² for each combination",
             x = "Ancestry", y = "Trait",
             fill = "Best Method") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave("prs_best_method_summary.pdf", p_best, width = 10, height = 6)

    # =========================================================================
    # 5. AUC COMPARISON (for binary traits)
    # =========================================================================
    if ("auc" %in% names(prs)) {
        prs_binary <- prs %>% filter(!is.na(auc))

        if (nrow(prs_binary) > 0) {
            p_auc <- ggplot(prs_binary, aes(x = method, y = auc, fill = ancestry)) +
                geom_bar(stat = "identity", position = position_dodge(0.8)) +
                geom_errorbar(aes(ymin = auc - auc_se, ymax = auc + auc_se),
                              position = position_dodge(0.8), width = 0.2) +
                geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
                facet_wrap(~trait, ncol = 2) +
                scale_fill_manual(values = ancestry_colors) +
                labs(title = "PRS AUC for Binary Traits",
                     x = "PRS Method", y = "AUC",
                     fill = "Ancestry") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))

            ggsave("prs_auc_comparison.pdf", p_auc, width = 12, height = 8)
        }
    }

    # =========================================================================
    # 6. PRS CORRELATION MATRIX (between methods)
    # =========================================================================
    # Calculate correlation between PRS methods within ancestry
    prs_corr <- prs %>%
        group_by(trait, ancestry) %>%
        summarize(
            n_methods = n(),
            mean_r2 = mean(r2),
            .groups = "drop"
        )

    p_corr <- ggplot(prs_corr, aes(x = ancestry, y = trait, fill = mean_r2, label = n_methods)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("n=%d\\nR²=%.3f", n_methods, mean_r2)),
                  color = "black", size = 2.5) +
        scale_fill_viridis_c() +
        labs(title = "PRS Methods Tested by Ancestry and Trait",
             subtitle = "Shows number of methods and mean R²",
             x = "Ancestry", y = "Trait") +
        theme_minimal()

    ggsave("prs_correlation_matrix.pdf", p_corr, width = 10, height = 6)

    # =========================================================================
    # 7. SUMMARY TABLE
    # =========================================================================
    summary_prs <- prs %>%
        group_by(trait, method) %>%
        summarize(
            n_ancestries = n(),
            mean_r2 = mean(r2, na.rm=TRUE),
            sd_r2 = sd(r2, na.rm=TRUE),
            max_r2 = max(r2, na.rm=TRUE),
            best_ancestry = ancestry[which.max(r2)],
            min_r2 = min(r2, na.rm=TRUE),
            worst_ancestry = ancestry[which.min(r2)],
            transferability_ratio = min_r2 / max_r2,  # Lower = worse transferability
            .groups = "drop"
        )

    write.table(summary_prs, "prs_summary_table.tsv", sep="\\t", row.names=FALSE, quote=FALSE)

    # Version info
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    ggplot2: "', packageVersion("ggplot2"), '"')
    ), "versions.yml")
    """
}

process CROSS_ANCESTRY_GWAS_PLOTS {
    tag "gwas_comparison"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2=3.4 conda-forge::r-dplyr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    path gwas_summaries      // Combined GWAS summaries by ancestry

    output:
    path "gwas_lambda_gc_comparison.pdf", emit: lambda
    path "gwas_n_significant_comparison.pdf", emit: n_sig
    path "gwas_effect_size_comparison.pdf", emit: effect_sizes
    path "gwas_shared_loci_heatmap.pdf", emit: shared_loci
    path "gwas_ancestry_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)

    # Read GWAS summaries
    gwas <- read.table("${gwas_summaries}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

    ancestry_colors <- c(
        "EUR" = "#4DAF4A", "AFR" = "#FF7F00", "EAS" = "#377EB8",
        "SAS" = "#984EA3", "AMR" = "#E41A1C", "MID" = "#A65628",
        "AAC" = "#F781BF", "AHI" = "#999999", "LATINO" = "#FFFF33"
    )

    # =========================================================================
    # 1. LAMBDA GC COMPARISON
    # =========================================================================
    if ("lambda_gc" %in% names(gwas)) {
        p_lambda <- ggplot(gwas, aes(x = ancestry, y = lambda_gc, fill = ancestry)) +
            geom_bar(stat = "identity", width = 0.7) +
            geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
            facet_wrap(~trait, ncol = 3) +
            scale_fill_manual(values = ancestry_colors) +
            labs(title = "Genomic Inflation Factor (λGC) by Ancestry",
                 subtitle = "Red line indicates no inflation (λ=1)",
                 x = "Ancestry", y = expression(lambda[GC])) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")

        ggsave("gwas_lambda_gc_comparison.pdf", p_lambda, width = 12, height = 8)
    }

    # =========================================================================
    # 2. NUMBER OF SIGNIFICANT HITS
    # =========================================================================
    if ("n_significant" %in% names(gwas)) {
        p_nsig <- ggplot(gwas, aes(x = ancestry, y = n_significant, fill = ancestry)) +
            geom_bar(stat = "identity", width = 0.7) +
            geom_text(aes(label = n_significant), vjust = -0.5, size = 3) +
            facet_wrap(~trait, scales = "free_y", ncol = 3) +
            scale_fill_manual(values = ancestry_colors) +
            labs(title = "Number of Genome-Wide Significant Variants by Ancestry",
                 subtitle = "P < 5×10⁻⁸",
                 x = "Ancestry", y = "N Significant SNPs") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")

        ggsave("gwas_n_significant_comparison.pdf", p_nsig, width = 12, height = 8)
    }

    # =========================================================================
    # 3. EFFECT SIZE COMPARISON (scatter plot of betas)
    # =========================================================================
    # This would require paired effect sizes - placeholder
    pdf("gwas_effect_size_comparison.pdf", width=10, height=8)
    plot.new()
    text(0.5, 0.5, "Effect size comparison requires paired variant data", cex=1.2)
    dev.off()

    # =========================================================================
    # 4. SHARED LOCI HEATMAP (Jaccard index)
    # =========================================================================
    pdf("gwas_shared_loci_heatmap.pdf", width=10, height=8)
    plot.new()
    text(0.5, 0.5, "Shared loci analysis requires individual variant data", cex=1.2)
    dev.off()

    # =========================================================================
    # 5. SUMMARY TABLE
    # =========================================================================
    summary_gwas <- gwas %>%
        group_by(trait) %>%
        summarize(
            n_ancestries = n(),
            total_significant = sum(n_significant, na.rm=TRUE),
            mean_lambda = mean(lambda_gc, na.rm=TRUE),
            max_significant_ancestry = ancestry[which.max(n_significant)],
            .groups = "drop"
        )

    write.table(summary_gwas, "gwas_ancestry_summary.tsv", sep="\\t", row.names=FALSE, quote=FALSE)

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

process SUMMARY_DASHBOARD {
    tag "dashboard"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-plotly conda-forge::r-htmlwidgets"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-plotly:4.10.3' :
        'quay.io/biocontainers/r-plotly:4.10.3' }"

    input:
    path gwas_summary
    path h2_summary
    path prs_summary
    path coloc_summary

    output:
    path "pipeline_summary_dashboard.html", emit: dashboard
    path "pipeline_summary_report.pdf", emit: pdf_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(plotly)
    library(htmlwidgets)
    library(dplyr)

    # Read all summaries
    gwas <- if(file.exists("${gwas_summary}")) read.table("${gwas_summary}", header=TRUE, sep="\\t") else data.frame()
    h2 <- if(file.exists("${h2_summary}")) read.table("${h2_summary}", header=TRUE, sep="\\t") else data.frame()
    prs <- if(file.exists("${prs_summary}")) read.table("${prs_summary}", header=TRUE, sep="\\t") else data.frame()
    coloc <- if(file.exists("${coloc_summary}")) read.table("${coloc_summary}", header=TRUE, sep="\\t") else data.frame()

    # Create interactive dashboard
    html_content <- '
    <!DOCTYPE html>
    <html>
    <head>
        <title>Ancestry-Aware GWAS Pipeline Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            .section { margin-bottom: 30px; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
            h1 { color: #2c3e50; }
            h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 5px; }
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #3498db; color: white; }
            tr:nth-child(even) { background-color: #f2f2f2; }
            .highlight { background-color: #e8f4f8; }
        </style>
    </head>
    <body>
        <h1>Ancestry-Aware GWAS Pipeline Summary Dashboard</h1>
        <p>Generated: ' + Sys.time() + '</p>

        <div class="section">
            <h2>GWAS Summary</h2>
            <p>Ancestry-stratified GWAS results</p>
        </div>

        <div class="section">
            <h2>Heritability Summary</h2>
            <p>SNP-heritability estimates by ancestry</p>
        </div>

        <div class="section">
            <h2>PRS Performance Summary</h2>
            <p>Polygenic risk score accuracy by method and ancestry</p>
        </div>

        <div class="section">
            <h2>Colocalization Summary</h2>
            <p>QTL colocalization results</p>
        </div>
    </body>
    </html>
    '

    writeLines(html_content, "pipeline_summary_dashboard.html")

    # Create PDF summary
    pdf("pipeline_summary_report.pdf", width=12, height=8)
    par(mfrow=c(2,2))
    plot.new()
    text(0.5, 0.5, "GWAS Summary", cex=2)
    plot.new()
    text(0.5, 0.5, "Heritability Summary", cex=2)
    plot.new()
    text(0.5, 0.5, "PRS Performance", cex=2)
    plot.new()
    text(0.5, 0.5, "Colocalization", cex=2)
    dev.off()

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}
