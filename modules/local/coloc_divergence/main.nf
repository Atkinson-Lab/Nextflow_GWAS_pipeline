/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SHARED vs DIVERGENT COLOCALIZATION ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Compares colocalization results between:
    - Ancestry-stratified GWAS (population-specific)
    - Meta-analysis GWAS (combined)

    Identifies:
    1. SHARED colocalizations: Present in meta-analysis AND multiple ancestries
       → Universal regulatory mechanisms, good drug targets
    2. DIVERGENT colocalizations: Present in only one ancestry
       → Population-specific regulatory effects, ancestry-aware medicine
    3. META-ONLY colocalizations: Present in meta but no single ancestry
       → May indicate heterogeneous effects that average out

    This helps identify population-specific vs universal drug targets.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COLOC_DIVERGENCE_ANALYSIS {
    tag "$trait"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-dplyr conda-forge::r-tidyr conda-forge::r-ggplot2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    path ancestry_coloc_files   // List of ancestry-stratified coloc results
    path meta_coloc_file        // Meta-analysis coloc result
    val trait                   // Trait name
    val pp4_threshold           // PP4 threshold for "colocalized" (default 0.8)

    output:
    path "${trait}.shared_coloc.tsv", emit: shared
    path "${trait}.divergent_coloc.tsv", emit: divergent
    path "${trait}.meta_only_coloc.tsv", emit: meta_only
    path "${trait}.coloc_comparison_summary.tsv", emit: summary
    path "${trait}.coloc_divergence.pdf", emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pp4_thresh = pp4_threshold ?: 0.8
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(tidyr)
    library(ggplot2)

    # =========================================================================
    # READ AND COMBINE ANCESTRY-STRATIFIED RESULTS
    # =========================================================================

    ancestry_files <- strsplit("${ancestry_coloc_files}", " ")[[1]]
    ancestry_files <- ancestry_files[file.exists(ancestry_files)]

    ancestry_results <- lapply(ancestry_files, function(f) {
        tryCatch({
            df <- read.table(f, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
            # Extract ancestry from filename or metadata
            ancestry <- gsub(".*\\\\.([A-Z]+)\\\\.coloc.*", "\\\\1", basename(f))
            if (ancestry == basename(f)) ancestry <- "UNKNOWN"
            df\$ancestry <- ancestry
            df
        }, error = function(e) NULL)
    })

    ancestry_df <- bind_rows(ancestry_results[!sapply(ancestry_results, is.null)])

    # =========================================================================
    # READ META-ANALYSIS RESULTS
    # =========================================================================

    meta_df <- tryCatch({
        read.table("${meta_coloc_file}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    }, error = function(e) {
        data.frame(gene = character(), qtl_type = character(), PP.H4.abf = numeric())
    })
    meta_df\$ancestry <- "META"

    # =========================================================================
    # DEFINE COLOCALIZED LOCI
    # =========================================================================

    pp4_threshold <- ${pp4_thresh}

    # Get colocalized gene-QTL pairs by ancestry
    ancestry_coloc <- ancestry_df %>%
        filter(PP.H4.abf >= pp4_threshold) %>%
        select(gene, qtl_type, ancestry, PP.H4.abf) %>%
        distinct()

    meta_coloc <- meta_df %>%
        filter(PP.H4.abf >= pp4_threshold) %>%
        select(gene, qtl_type, PP.H4.abf) %>%
        mutate(ancestry = "META") %>%
        distinct()

    # =========================================================================
    # CLASSIFY: SHARED vs DIVERGENT vs META-ONLY
    # =========================================================================

    # Create gene-QTL identifiers
    ancestry_coloc <- ancestry_coloc %>%
        mutate(locus = paste(gene, qtl_type, sep = ":"))

    meta_coloc <- meta_coloc %>%
        mutate(locus = paste(gene, qtl_type, sep = ":"))

    # Count how many ancestries show colocalization per locus
    ancestry_counts <- ancestry_coloc %>%
        group_by(locus) %>%
        summarize(
            n_ancestries = n_distinct(ancestry),
            ancestries = paste(sort(unique(ancestry)), collapse = ","),
            mean_pp4 = mean(PP.H4.abf),
            max_pp4 = max(PP.H4.abf),
            .groups = "drop"
        )

    # SHARED: In meta-analysis AND >= 2 ancestries
    shared_loci <- ancestry_counts %>%
        filter(n_ancestries >= 2, locus %in% meta_coloc\$locus) %>%
        mutate(category = "SHARED")

    # DIVERGENT: Only in 1 ancestry (population-specific)
    divergent_loci <- ancestry_counts %>%
        filter(n_ancestries == 1) %>%
        mutate(category = "DIVERGENT")

    # META-ONLY: In meta but not reaching threshold in any single ancestry
    meta_only_loci <- meta_coloc %>%
        filter(!locus %in% ancestry_coloc\$locus) %>%
        select(locus, PP.H4.abf) %>%
        rename(meta_pp4 = PP.H4.abf) %>%
        mutate(category = "META_ONLY")

    # =========================================================================
    # OUTPUT RESULTS
    # =========================================================================

    # Shared colocalizations (universal mechanisms)
    shared_out <- shared_loci %>%
        separate(locus, into = c("gene", "qtl_type"), sep = ":", remove = FALSE) %>%
        left_join(meta_coloc %>% select(locus, PP.H4.abf) %>% rename(meta_pp4 = PP.H4.abf), by = "locus")

    write.table(shared_out, "${trait}.shared_coloc.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    # Divergent colocalizations (population-specific)
    divergent_out <- divergent_loci %>%
        separate(locus, into = c("gene", "qtl_type"), sep = ":", remove = FALSE)

    write.table(divergent_out, "${trait}.divergent_coloc.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    # Meta-only colocalizations
    meta_only_out <- meta_only_loci %>%
        separate(locus, into = c("gene", "qtl_type"), sep = ":", remove = FALSE)

    write.table(meta_only_out, "${trait}.meta_only_coloc.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    # Summary statistics
    summary_df <- data.frame(
        trait = "${trait}",
        n_shared = nrow(shared_loci),
        n_divergent = nrow(divergent_loci),
        n_meta_only = nrow(meta_only_loci),
        n_total_coloc = nrow(ancestry_counts),
        pct_shared = round(100 * nrow(shared_loci) / max(1, nrow(ancestry_counts)), 1),
        pct_divergent = round(100 * nrow(divergent_loci) / max(1, nrow(ancestry_counts)), 1),
        divergent_ancestries = paste(unique(divergent_loci\$ancestries), collapse = "; ")
    )

    write.table(summary_df, "${trait}.coloc_comparison_summary.tsv",
                sep = "\\t", row.names = FALSE, quote = FALSE)

    # =========================================================================
    # VISUALIZATION
    # =========================================================================

    # Combine for plotting
    all_loci <- bind_rows(
        shared_out %>% mutate(category = "Shared"),
        divergent_out %>% mutate(category = paste0("Divergent (", ancestries, ")")),
        meta_only_out %>% mutate(category = "Meta-only")
    )

    if (nrow(all_loci) > 0) {
        # Bar plot of categories
        p1 <- ggplot(all_loci, aes(x = category, fill = category)) +
            geom_bar() +
            scale_fill_manual(values = c(
                "Shared" = "#2ecc71",
                "Meta-only" = "#3498db",
                "Divergent" = "#e74c3c"
            ), na.value = "#e74c3c") +
            labs(title = paste("Colocalization Categories:", "${trait}"),
                 subtitle = paste("PP4 threshold:", pp4_threshold),
                 x = "Category",
                 y = "Number of Gene-QTL Pairs") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")

        # If we have ancestry-specific data, show PP4 by ancestry
        if (nrow(ancestry_df) > 0 && "PP.H4.abf" %in% colnames(ancestry_df)) {
            ancestry_summary <- ancestry_df %>%
                filter(!is.na(PP.H4.abf)) %>%
                group_by(ancestry) %>%
                summarize(
                    n_tested = n(),
                    n_coloc = sum(PP.H4.abf >= pp4_threshold),
                    mean_pp4 = mean(PP.H4.abf),
                    .groups = "drop"
                )

            p2 <- ggplot(ancestry_summary, aes(x = reorder(ancestry, -n_coloc), y = n_coloc, fill = ancestry)) +
                geom_bar(stat = "identity") +
                labs(title = "Colocalizations by Ancestry",
                     x = "Ancestry",
                     y = paste("N with PP4 >=", pp4_threshold)) +
                theme_bw() +
                theme(legend.position = "none")
        } else {
            p2 <- ggplot() + theme_void() +
                annotate("text", x = 0.5, y = 0.5, label = "No ancestry data")
        }

        # Combine plots
        library(patchwork)
        combined <- p1 + p2 + plot_layout(ncol = 2)

        ggsave("${trait}.coloc_divergence.pdf", combined, width = 12, height = 6)
    } else {
        pdf("${trait}.coloc_divergence.pdf", width = 8, height = 6)
        plot.new()
        text(0.5, 0.5, "No colocalizations found", cex = 2)
        dev.off()
    }

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

process COLOC_ANCESTRY_HEATMAP {
    tag "$trait"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-pheatmap conda-forge::r-dplyr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    path coloc_files    // All coloc results (ancestry + meta)
    val trait

    output:
    path "${trait}.coloc_ancestry_heatmap.pdf", emit: heatmap
    path "${trait}.coloc_pp4_matrix.tsv", emit: matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(tidyr)
    library(pheatmap)

    # Read all coloc files
    coloc_files <- strsplit("${coloc_files}", " ")[[1]]
    coloc_files <- coloc_files[file.exists(coloc_files)]

    all_results <- lapply(coloc_files, function(f) {
        tryCatch({
            df <- read.table(f, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
            # Extract ancestry/analysis type from filename
            if (grepl("meta", tolower(basename(f)))) {
                analysis <- "META"
            } else {
                analysis <- gsub(".*\\\\.([A-Z]+)\\\\.coloc.*", "\\\\1", basename(f))
                if (analysis == basename(f)) analysis <- "UNKNOWN"
            }
            df\$analysis <- analysis
            df
        }, error = function(e) NULL)
    })

    combined <- bind_rows(all_results[!sapply(all_results, is.null)])

    if (nrow(combined) > 0 && "PP.H4.abf" %in% colnames(combined)) {
        # Create gene-QTL x ancestry matrix of PP4 values
        pp4_wide <- combined %>%
            mutate(locus = paste(gene, qtl_type, sep = ":")) %>%
            select(locus, analysis, PP.H4.abf) %>%
            group_by(locus, analysis) %>%
            summarize(PP4 = max(PP.H4.abf), .groups = "drop") %>%
            pivot_wider(names_from = analysis, values_from = PP4, values_fill = 0) %>%
            as.data.frame()

        rownames(pp4_wide) <- pp4_wide\$locus
        pp4_wide\$locus <- NULL

        # Save matrix
        write.table(pp4_wide, "${trait}.coloc_pp4_matrix.tsv",
                    sep = "\\t", row.names = TRUE, quote = FALSE)

        # Create heatmap
        # Limit to top loci if too many
        if (nrow(pp4_wide) > 100) {
            top_loci <- names(sort(rowSums(pp4_wide), decreasing = TRUE)[1:100])
            pp4_wide <- pp4_wide[top_loci, , drop = FALSE]
        }

        pdf("${trait}.coloc_ancestry_heatmap.pdf", width = 10, height = max(8, nrow(pp4_wide) * 0.15))
        pheatmap(as.matrix(pp4_wide),
                 main = paste("Colocalization PP4 by Ancestry:", "${trait}"),
                 color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
                 breaks = seq(0, 1, length.out = 101),
                 cluster_rows = nrow(pp4_wide) > 2,
                 cluster_cols = ncol(pp4_wide) > 2,
                 show_rownames = nrow(pp4_wide) <= 50,
                 fontsize_row = 8,
                 annotation_legend = TRUE)
        dev.off()
    } else {
        # Empty outputs
        write.table(data.frame(), "${trait}.coloc_pp4_matrix.tsv", sep = "\\t")
        pdf("${trait}.coloc_ancestry_heatmap.pdf", width = 8, height = 6)
        plot.new()
        text(0.5, 0.5, "No colocalization data", cex = 2)
        dev.off()
    }

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTERPRETATION GUIDE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    SHARED COLOCALIZATIONS:
    - Gene-QTL pairs with PP4 >= threshold in META-analysis AND >= 2 ancestries
    - Represent UNIVERSAL regulatory mechanisms
    - High-confidence targets for drug development
    - Effects are consistent across populations

    DIVERGENT COLOCALIZATIONS:
    - Gene-QTL pairs colocalized in ONLY ONE ancestry
    - Represent POPULATION-SPECIFIC regulatory effects
    - Important for precision medicine / ancestry-aware treatment
    - May reflect:
      a) Different LD patterns across populations
      b) Different allele frequencies
      c) True biological differences in gene regulation
      d) Population-specific environmental interactions

    META-ONLY COLOCALIZATIONS:
    - Gene-QTL pairs colocalized in META but not any single ancestry
    - May indicate:
      a) Weak effects that only reach significance with combined power
      b) Heterogeneous effects across populations
      c) Statistical artifacts
    - Require careful interpretation

    CLINICAL IMPLICATIONS:
    - SHARED: Safe to pursue as universal targets
    - DIVERGENT: Consider for ancestry-specific therapies
    - META-ONLY: Investigate further before clinical translation

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
