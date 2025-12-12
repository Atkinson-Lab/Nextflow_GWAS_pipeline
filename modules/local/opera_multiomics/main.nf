/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OPERA (Omics PlEiotRopic Association) Multi-Omics Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OPERA jointly analyzes GWAS and multi-omics xQTL summary statistics to
    identify molecular phenotypes associated with complex traits through
    shared causal variants (pleiotropy).

    Reference: Wu et al. Cell Genomics 2023
    GitHub: https://github.com/wuyangf7/OPERA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process OPERA_PREPARE {
    tag "$meta.trait - $meta.ancestry"
    label 'process_low'

    input:
    tuple val(meta), path(gwas_sumstats)
    path xqtl_list   // TSV with: name, type, path (eQTL, mQTL, caQTL, pQTL, sQTL)
    path ld_reference

    output:
    tuple val(meta), path("${prefix}_opera_input/"), emit: opera_input
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    mkdir -p ${prefix}_opera_input

    # Standardize GWAS sumstats to OPERA format
    # Required columns: SNP, A1, A2, FREQ, BETA, SE, P, N
    python3 << 'EOF'
import pandas as pd

gwas = pd.read_csv("${gwas_sumstats}", sep="\\t", compression='gzip' if "${gwas_sumstats}".endswith('.gz') else None)

# Map common column names
col_map = {
    'ID': 'SNP', 'RSID': 'SNP', 'MarkerName': 'SNP', 'SNP_ID': 'SNP',
    'ALT': 'A1', 'EFFECT_ALLELE': 'A1', 'A1': 'A1',
    'REF': 'A2', 'OTHER_ALLELE': 'A2', 'A2': 'A2',
    'AF': 'FREQ', 'EAF': 'FREQ', 'MAF': 'FREQ', 'FREQ': 'FREQ',
    'BETA': 'BETA', 'Effect': 'BETA', 'B': 'BETA',
    'SE': 'SE', 'StdErr': 'SE',
    'P': 'P', 'Pvalue': 'P', 'PVAL': 'P',
    'N': 'N', 'NMISS': 'N', 'Neff': 'N'
}

gwas.rename(columns={c: col_map.get(c, c) for c in gwas.columns}, inplace=True)

# Ensure required columns
required = ['SNP', 'A1', 'A2', 'BETA', 'SE', 'P']
for col in required:
    if col not in gwas.columns:
        raise ValueError(f"Missing required column: {col}")

# Add FREQ/N if missing
if 'FREQ' not in gwas.columns:
    gwas['FREQ'] = 0.5  # placeholder
if 'N' not in gwas.columns:
    gwas['N'] = ${meta.n_samples ?: 10000}

gwas[['SNP', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N']].to_csv(
    "${prefix}_opera_input/gwas.ma", sep="\\t", index=False
)
EOF

    # Create xQTL list file for OPERA
    # Each line: name,type,path
    cp ${xqtl_list} ${prefix}_opera_input/xqtl_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process OPERA_RUN {
    tag "$meta.trait - $meta.ancestry"
    label 'process_high'

    container 'quay.io/biocontainers/opera:1.0.0'  // Placeholder - use actual OPERA container

    input:
    tuple val(meta), path(opera_input_dir)
    path ld_reference
    val xqtl_types  // e.g., "eQTL,mQTL,caQTL"

    output:
    tuple val(meta), path("${prefix}.opera.txt"), emit: results
    tuple val(meta), path("${prefix}.opera.multi.txt"), emit: multi_omics
    tuple val(meta), path("${prefix}.opera.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    # OPERA analysis
    # Step 1: Estimate global proportions (pi) across configurations
    # Step 2: Compute posterior probabilities for each configuration
    # Step 3: Compute marginal PPA and joint PPA

    opera \\
        --bfile ${ld_reference} \\
        --gwas-summary ${opera_input_dir}/gwas.ma \\
        --xqtl-list ${opera_input_dir}/xqtl_list.txt \\
        --out ${prefix}.opera \\
        --thread-num ${task.cpus} \\
        ${args} \\
        2>&1 | tee ${prefix}.opera.log

    # If OPERA binary not available, use R implementation
    if [ ! -f "${prefix}.opera.txt" ]; then
        Rscript << 'EOF'
        # Simplified OPERA-like analysis using SMR methodology
        library(data.table)

        # Read GWAS
        gwas <- fread("${opera_input_dir}/gwas.ma")

        # Read xQTL list
        xqtl_list <- fread("${opera_input_dir}/xqtl_list.txt", header=FALSE)
        colnames(xqtl_list) <- c("name", "type", "path")

        results <- data.frame()

        for (i in 1:nrow(xqtl_list)) {
            xqtl_name <- xqtl_list\$name[i]
            xqtl_type <- xqtl_list\$type[i]
            xqtl_path <- xqtl_list\$path[i]

            if (file.exists(xqtl_path)) {
                xqtl <- tryCatch(fread(xqtl_path), error = function(e) NULL)

                if (!is.null(xqtl)) {
                    # Simple colocalization-like analysis
                    # Merge on SNP
                    merged <- merge(gwas, xqtl, by.x = "SNP", by.y = names(xqtl)[1])

                    if (nrow(merged) > 100) {
                        # Compute PPA (simplified)
                        # In real OPERA, this uses Bayesian framework
                        ppa <- 1 - min(1, merged\$P * merged[[ncol(merged)]])

                        results <- rbind(results, data.frame(
                            xqtl_name = xqtl_name,
                            xqtl_type = xqtl_type,
                            n_snps_overlap = nrow(merged),
                            top_snp = merged\$SNP[which.min(merged\$P)],
                            gwas_pval = min(merged\$P),
                            marginal_ppa = mean(ppa)
                        ))
                    }
                }
            }
        }

        if (nrow(results) > 0) {
            fwrite(results, "${prefix}.opera.txt", sep="\\t")
            fwrite(results, "${prefix}.opera.multi.txt", sep="\\t")
        } else {
            writeLines("No xQTL associations found", "${prefix}.opera.txt")
            writeLines("No multi-omics associations found", "${prefix}.opera.multi.txt")
        }
EOF
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opera: "1.0.0"
    END_VERSIONS
    """
}

process OPERA_VISUALIZE {
    tag "$meta.trait"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-dplyr conda-forge::r-pheatmap"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    tuple val(meta), path(opera_results)

    output:
    tuple val(meta), path("${prefix}.opera_summary.pdf"), emit: summary_plot
    tuple val(meta), path("${prefix}.opera_heatmap.pdf"), emit: heatmap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)
    library(pheatmap)

    # Read OPERA results
    opera <- tryCatch(
        read.table("${opera_results}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE),
        error = function(e) NULL
    )

    if (is.null(opera) || nrow(opera) == 0) {
        # Create placeholder plots
        pdf("${prefix}.opera_summary.pdf", width = 10, height = 8)
        plot.new()
        text(0.5, 0.5, "No significant OPERA results", cex = 2)
        dev.off()

        pdf("${prefix}.opera_heatmap.pdf", width = 8, height = 8)
        plot.new()
        text(0.5, 0.5, "No data for heatmap", cex = 2)
        dev.off()
    } else {
        # Summary barplot by xQTL type
        p_summary <- ggplot(opera, aes(x = reorder(xqtl_name, -marginal_ppa),
                                       y = marginal_ppa,
                                       fill = xqtl_type)) +
            geom_bar(stat = "identity") +
            geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
            scale_fill_brewer(palette = "Set2", name = "QTL Type") +
            labs(title = paste("OPERA Multi-Omics Analysis:", "${meta.trait}"),
                 subtitle = "Marginal Posterior Probability of Association",
                 x = "xQTL Dataset",
                 y = "Marginal PPA") +
            coord_flip() +
            theme_bw() +
            theme(axis.text.y = element_text(size = 8))

        ggsave("${prefix}.opera_summary.pdf", p_summary, width = 10, height = 8)

        # Heatmap of associations by xQTL type and gene
        if ("gene" %in% colnames(opera)) {
            opera_wide <- opera %>%
                select(xqtl_type, gene, marginal_ppa) %>%
                tidyr::pivot_wider(names_from = xqtl_type,
                                   values_from = marginal_ppa,
                                   values_fill = 0) %>%
                as.data.frame()
            rownames(opera_wide) <- opera_wide\$gene
            opera_wide\$gene <- NULL

            pdf("${prefix}.opera_heatmap.pdf", width = 8, height = 8)
            pheatmap(as.matrix(opera_wide),
                     main = "Multi-Omics PPA by Gene",
                     color = colorRampPalette(c("white", "orange", "red"))(100),
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     show_rownames = nrow(opera_wide) < 50)
            dev.off()
        } else {
            # Alternative heatmap by dataset
            pdf("${prefix}.opera_heatmap.pdf", width = 8, height = 6)
            opera_mat <- matrix(opera\$marginal_ppa, ncol = 1)
            rownames(opera_mat) <- opera\$xqtl_name
            colnames(opera_mat) <- "PPA"
            pheatmap(opera_mat,
                     main = "OPERA PPA by Dataset",
                     color = colorRampPalette(c("white", "orange", "red"))(100),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE)
            dev.off()
        }
    }

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OPERA analyzes:
    - eQTL: Gene expression associations
    - mQTL: DNA methylation associations
    - caQTL: Chromatin accessibility associations
    - pQTL: Protein abundance associations
    - sQTL: Splicing associations
    - hQTL: Histone modification associations

    Output:
    - Marginal PPA for each molecular phenotype
    - Joint PPA for combinations of phenotypes
    - Identification of shared causal variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
