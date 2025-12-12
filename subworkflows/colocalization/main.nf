/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLOCALIZATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Colocalization analysis with eQTL/sQTL/pQTL/mQTL/caQTL/hQTL datasets

    IMPORTANT: This workflow runs colocalization on BOTH:
    1. Ancestry-stratified GWAS results (per-ancestry analysis)
    2. Meta-analysis results (combined across ancestries)

    The rationale:
    - Ancestry-stratified: Identifies population-specific regulatory mechanisms
    - Meta-analysis: Higher power to detect shared causal variants

    Supports pooling of diverse QTL datasets for improved representation
    across ancestries (addressing GTEx's limited diversity).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_QTL_DATA           } from '../../modules/local/prepare_qtl_data'
include { POOL_QTL_DATASETS          } from '../../modules/local/pool_qtl_datasets'
include { COLOC                      } from '../../modules/local/coloc'
include { ECAVIAR                    } from '../../modules/local/ecaviar'
include { FASTENLOC                  } from '../../modules/local/fastenloc'
include { COLOC_SUMMARY              } from '../../modules/local/coloc_summary'
include { COLOC_DIVERGENCE_ANALYSIS  } from '../../modules/local/coloc_divergence'
include { COLOC_ANCESTRY_HEATMAP     } from '../../modules/local/coloc_divergence'

workflow COLOCALIZATION {
    take:
    ch_ancestry_gwas   // channel: [ meta, sumstats ] - Ancestry-stratified GWAS results
    ch_meta_gwas       // channel: [ meta, sumstats ] - Meta-analysis GWAS results
    ch_qtl_data        // channel: [ name, type, qtl_file ]
    coloc_method       // value: 'coloc', 'ecaviar', 'fastenloc'
    p1                 // value: prior for trait association
    p2                 // value: prior for QTL association
    p12                // value: prior for colocalization

    main:
    ch_versions = Channel.empty()
    ch_coloc_ancestry = Channel.empty()
    ch_coloc_meta = Channel.empty()

    // =========================================================================
    // PREPARE QTL DATA
    // =========================================================================

    // Standardize QTL data format
    PREPARE_QTL_DATA(
        ch_qtl_data
    )
    ch_versions = ch_versions.mix(PREPARE_QTL_DATA.out.versions)

    // Optionally pool diverse QTL datasets to create mega-set
    // This addresses the limited diversity in GTEx by combining multiple sources
    ch_qtl_by_type = PREPARE_QTL_DATA.out.standardized
        .map { name, type, qtl_file -> [type, name, qtl_file] }
        .groupTuple(by: 0)

    POOL_QTL_DATASETS(
        ch_qtl_by_type
    )
    ch_versions = ch_versions.mix(POOL_QTL_DATASETS.out.versions)

    // Use pooled datasets for colocalization
    ch_qtl_pooled = POOL_QTL_DATASETS.out.pooled

    // =========================================================================
    // COLOCALIZATION: ANCESTRY-STRATIFIED GWAS
    // =========================================================================
    // Run colocalization on each ancestry-specific GWAS separately
    // This identifies population-specific regulatory mechanisms
    // Best practice: Match QTL ancestry to GWAS ancestry when possible

    // Tag ancestry results
    ch_ancestry_gwas_tagged = ch_ancestry_gwas
        .map { meta, sumstats ->
            def new_meta = meta.clone()
            new_meta.analysis_type = 'ancestry_stratified'
            [new_meta, sumstats]
        }

    // Create GWAS x QTL combinations for ancestry-stratified analysis
    ch_ancestry_coloc_input = ch_ancestry_gwas_tagged.combine(ch_qtl_pooled)

    // =========================================================================
    // COLOCALIZATION: META-ANALYSIS GWAS
    // =========================================================================
    // Run colocalization on meta-analysis results
    // Higher power due to larger sample size
    // Identifies shared causal variants across ancestries

    // Tag meta-analysis results
    ch_meta_gwas_tagged = ch_meta_gwas
        .map { meta, sumstats ->
            def new_meta = meta.clone()
            new_meta.analysis_type = 'meta_analysis'
            [new_meta, sumstats]
        }

    // Create GWAS x QTL combinations for meta-analysis
    ch_meta_coloc_input = ch_meta_gwas_tagged.combine(ch_qtl_pooled)

    // =========================================================================
    // COMBINE AND RUN COLOCALIZATION
    // =========================================================================
    // Merge both ancestry-stratified and meta-analysis inputs

    ch_all_coloc_input = ch_ancestry_coloc_input.mix(ch_meta_coloc_input)

    if (coloc_method == 'coloc') {
        COLOC(
            ch_all_coloc_input,
            p1,
            p2,
            p12
        )
        ch_coloc_all = COLOC.out.results
        ch_versions = ch_versions.mix(COLOC.out.versions)

    } else if (coloc_method == 'ecaviar') {
        ECAVIAR(
            ch_all_coloc_input
        )
        ch_coloc_all = ECAVIAR.out.results
        ch_versions = ch_versions.mix(ECAVIAR.out.versions)

    } else if (coloc_method == 'fastenloc') {
        FASTENLOC(
            ch_all_coloc_input
        )
        ch_coloc_all = FASTENLOC.out.results
        ch_versions = ch_versions.mix(FASTENLOC.out.versions)
    }

    // =========================================================================
    // SPLIT RESULTS BY ANALYSIS TYPE
    // =========================================================================

    // Separate ancestry-stratified and meta-analysis results
    ch_coloc_ancestry = ch_coloc_all
        .filter { meta, results -> meta.analysis_type == 'ancestry_stratified' }

    ch_coloc_meta = ch_coloc_all
        .filter { meta, results -> meta.analysis_type == 'meta_analysis' }

    // =========================================================================
    // SUMMARIZE COLOCALIZATION RESULTS
    // =========================================================================
    // Generate comprehensive summary across all analyses

    COLOC_SUMMARY(
        ch_coloc_all.collect { it[1] }
    )
    ch_versions = ch_versions.mix(COLOC_SUMMARY.out.versions)

    // =========================================================================
    // SHARED vs DIVERGENT COLOCALIZATION ANALYSIS
    // =========================================================================
    // Compare ancestry-stratified to meta-analysis to identify:
    // - SHARED: Universal regulatory mechanisms (in meta + multiple ancestries)
    // - DIVERGENT: Population-specific effects (only in one ancestry)
    // - META-ONLY: Effects only significant when combined

    // Group coloc results by trait for divergence analysis
    ch_ancestry_by_trait = ch_coloc_ancestry
        .map { meta, results -> [meta.trait, results] }
        .groupTuple()

    ch_meta_by_trait = ch_coloc_meta
        .map { meta, results -> [meta.trait, results] }

    // Join ancestry and meta results by trait
    ch_divergence_input = ch_ancestry_by_trait
        .join(ch_meta_by_trait)
        .map { trait, ancestry_files, meta_file ->
            [ancestry_files, meta_file, trait]
        }

    COLOC_DIVERGENCE_ANALYSIS(
        ch_divergence_input.map { it[0] },  // ancestry files
        ch_divergence_input.map { it[1] },  // meta file
        ch_divergence_input.map { it[2] },  // trait name
        0.8  // PP4 threshold
    )
    ch_versions = ch_versions.mix(COLOC_DIVERGENCE_ANALYSIS.out.versions)

    // Heatmap of PP4 across ancestries
    ch_all_coloc_files = ch_coloc_all
        .map { meta, results -> [meta.trait, results] }
        .groupTuple()

    COLOC_ANCESTRY_HEATMAP(
        ch_all_coloc_files.map { trait, files -> files },
        ch_all_coloc_files.map { trait, files -> trait }
    )
    ch_versions = ch_versions.mix(COLOC_ANCESTRY_HEATMAP.out.versions)

    emit:
    // All colocalization results combined
    coloc_results         = ch_coloc_all               // channel: [ meta, coloc_results ]
    // Ancestry-stratified results only
    coloc_ancestry        = ch_coloc_ancestry          // channel: [ meta, coloc_results ]
    // Meta-analysis results only
    coloc_meta            = ch_coloc_meta              // channel: [ meta, coloc_results ]
    // Shared vs Divergent analysis
    shared_coloc          = COLOC_DIVERGENCE_ANALYSIS.out.shared     // channel: [ shared_coloc.tsv ]
    divergent_coloc       = COLOC_DIVERGENCE_ANALYSIS.out.divergent  // channel: [ divergent_coloc.tsv ]
    meta_only_coloc       = COLOC_DIVERGENCE_ANALYSIS.out.meta_only  // channel: [ meta_only_coloc.tsv ]
    divergence_summary    = COLOC_DIVERGENCE_ANALYSIS.out.summary    // channel: [ summary.tsv ]
    divergence_plot       = COLOC_DIVERGENCE_ANALYSIS.out.plot       // channel: [ divergence.pdf ]
    ancestry_heatmap      = COLOC_ANCESTRY_HEATMAP.out.heatmap       // channel: [ heatmap.pdf ]
    // Summary statistics
    summary               = COLOC_SUMMARY.out.summary  // channel: [ summary_file ]
    pooled_qtl            = ch_qtl_pooled              // channel: [ type, pooled_qtl_file ]
    versions              = ch_versions                // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USAGE NOTES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    To use this workflow in main.nf:

    COLOCALIZATION(
        GWAS.out.ancestry_results,    // Per-ancestry GWAS (e.g., EUR, AFR, EAS separately)
        META_ANALYSIS.out.results,    // Combined meta-analysis results
        ch_qtl_datasets,              // QTL data from curated sources
        params.coloc_method,
        params.coloc_p1,
        params.coloc_p2,
        params.coloc_p12
    )

    Results interpretation:
    - coloc_ancestry: Population-specific regulatory effects
      Use for identifying ancestry-specific drug targets
    - coloc_meta: Shared effects across populations
      Higher confidence for universal mechanisms

    For best results with diverse cohorts:
    - Match QTL datasets to GWAS ancestry when possible
    - Use multi-ancestry QTL sources (MESA, PAGE, etc.)
    - Include ancestry-specific QTLs (GENOA for AFR, BBJ for EAS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
