/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VISUALIZATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generate publication-quality plots
    Manhattan, QQ, Regional association, Miami, Forest plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MANHATTAN_PLOT      } from '../../modules/local/manhattan_plot'
include { QQ_PLOT             } from '../../modules/local/qq_plot'
include { REGIONAL_PLOT       } from '../../modules/local/regional_plot'
include { MIAMI_PLOT          } from '../../modules/local/miami_plot'
include { FOREST_PLOT         } from '../../modules/local/forest_plot'
include { PRS_PLOT            } from '../../modules/local/prs_plot'
include { LAMBDA_GC           } from '../../modules/local/lambda_gc'

workflow VISUALIZATION {
    take:
    ch_gwas_results       // channel: [ meta, sumstats ]
    ch_significant        // channel: [ meta, sig_variants ]
    gene_annotation       // file: gene annotation file

    main:
    ch_versions = Channel.empty()
    ch_plots = Channel.empty()

    // =========================================================================
    // CALCULATE GENOMIC INFLATION FACTOR
    // =========================================================================

    LAMBDA_GC(
        ch_gwas_results
    )
    ch_versions = ch_versions.mix(LAMBDA_GC.out.versions)

    // =========================================================================
    // MANHATTAN PLOTS
    // =========================================================================

    // Per-ancestry Manhattan plots
    MANHATTAN_PLOT(
        ch_gwas_results.join(LAMBDA_GC.out.lambda)
    )
    ch_plots = ch_plots.mix(MANHATTAN_PLOT.out.plots)
    ch_versions = ch_versions.mix(MANHATTAN_PLOT.out.versions)

    // =========================================================================
    // QQ PLOTS
    // =========================================================================

    QQ_PLOT(
        ch_gwas_results.join(LAMBDA_GC.out.lambda)
    )
    ch_plots = ch_plots.mix(QQ_PLOT.out.plots)
    ch_versions = ch_versions.mix(QQ_PLOT.out.versions)

    // =========================================================================
    // REGIONAL ASSOCIATION PLOTS (LocusZoom-style)
    // =========================================================================

    // Extract top loci for regional plots
    ch_top_loci = ch_significant
        .flatMap { meta, sig_file ->
            // Return channel of [meta, locus_info] for each significant locus
            // This would be parsed from the sig_file
            [[meta, sig_file]]
        }

    REGIONAL_PLOT(
        ch_top_loci,
        gene_annotation
    )
    ch_plots = ch_plots.mix(REGIONAL_PLOT.out.plots)
    ch_versions = ch_versions.mix(REGIONAL_PLOT.out.versions)

    // =========================================================================
    // MIAMI PLOTS (comparing ancestries)
    // =========================================================================

    // Create ancestry pairs for Miami plots
    ch_for_miami = ch_gwas_results
        .map { meta, sumstats -> [[trait: meta.trait], meta.ancestry, sumstats] }
        .groupTuple(by: 0)
        .filter { trait_meta, ancestries, sumstats_list ->
            ancestries.size() >= 2
        }
        .flatMap { trait_meta, ancestries, sumstats_list ->
            def pairs = []
            // Create pairs comparing each ancestry to EUR (if available) or first ancestry
            def ref_idx = ancestries.indexOf('EUR') >= 0 ? ancestries.indexOf('EUR') : 0
            for (int i = 0; i < ancestries.size(); i++) {
                if (i != ref_idx) {
                    pairs << [
                        trait_meta + [anc1: ancestries[ref_idx], anc2: ancestries[i]],
                        sumstats_list[ref_idx],
                        sumstats_list[i]
                    ]
                }
            }
            return pairs
        }

    MIAMI_PLOT(
        ch_for_miami
    )
    ch_plots = ch_plots.mix(MIAMI_PLOT.out.plots)
    ch_versions = ch_versions.mix(MIAMI_PLOT.out.versions)

    emit:
    plots        = ch_plots                    // channel: [ plot_files ]
    lambda_gc    = LAMBDA_GC.out.lambda        // channel: [ meta, lambda ]
    versions     = ch_versions                 // channel: [ versions.yml ]
}
