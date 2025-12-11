/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HERITABILITY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Heritability estimation overall and by ancestry
    Includes GxE heritability assessment and cross-ancestry genetic correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LDSC_MUNGE        } from '../../modules/local/ldsc_munge'
include { LDSC_H2           } from '../../modules/local/ldsc_h2'
include { LDSC_RG           } from '../../modules/local/ldsc_rg'
include { GCTA_GREML        } from '../../modules/local/gcta_greml'
include { GCTA_GXE          } from '../../modules/local/gcta_gxe'
include { BOLT_REML         } from '../../modules/local/bolt_reml'
include { H2_COMPARISON     } from '../../modules/local/h2_comparison'
include { GXE_ANALYSIS      } from '../../modules/local/gxe_analysis'

workflow HERITABILITY {
    take:
    ch_gwas_results        // channel: [ meta, sumstats ]
    ch_meta_results        // channel: [ meta, meta_sumstats ]
    h2_method              // value: 'ldsc', 'gcta', 'bolt-reml'
    ancestry_specific_h2   // boolean
    gxe_heritability       // boolean
    genetic_correlation    // boolean
    ld_reference_dir       // file: LD reference directory

    main:
    ch_versions = Channel.empty()
    ch_h2_results = Channel.empty()
    ch_rg_results = Channel.empty()
    ch_gxe_results = Channel.empty()

    // =========================================================================
    // ANCESTRY-SPECIFIC HERITABILITY
    // =========================================================================

    if (ancestry_specific_h2) {
        if (h2_method == 'ldsc') {
            // Munge sumstats for LDSC
            LDSC_MUNGE(
                ch_gwas_results
            )
            ch_versions = ch_versions.mix(LDSC_MUNGE.out.versions)

            // Calculate h² per ancestry
            LDSC_H2(
                LDSC_MUNGE.out.munged,
                ld_reference_dir
            )
            ch_h2_results = ch_h2_results.mix(LDSC_H2.out.h2)
            ch_versions = ch_versions.mix(LDSC_H2.out.versions)

        } else if (h2_method == 'gcta') {
            GCTA_GREML(
                ch_gwas_results
            )
            ch_h2_results = ch_h2_results.mix(GCTA_GREML.out.h2)
            ch_versions = ch_versions.mix(GCTA_GREML.out.versions)

        } else if (h2_method == 'bolt-reml') {
            BOLT_REML(
                ch_gwas_results
            )
            ch_h2_results = ch_h2_results.mix(BOLT_REML.out.h2)
            ch_versions = ch_versions.mix(BOLT_REML.out.versions)
        }
    }

    // =========================================================================
    // OVERALL HERITABILITY (from meta-analysis)
    // =========================================================================

    if (ch_meta_results && h2_method == 'ldsc') {
        LDSC_MUNGE.mix(ch_meta_results).set { ch_meta_to_munge }

        // Calculate overall h² from meta-analysis results
        ch_meta_munged = ch_meta_results
            .map { meta, sumstats -> [meta + [ancestry: 'META'], sumstats] }

        LDSC_MUNGE(ch_meta_munged)

        LDSC_H2(
            LDSC_MUNGE.out.munged,
            ld_reference_dir
        )
        ch_h2_results = ch_h2_results.mix(LDSC_H2.out.h2)
    }

    // =========================================================================
    // GxE HERITABILITY COMPONENT
    // =========================================================================

    if (gxe_heritability) {
        // Estimate GxE heritability to assess environmental interaction
        // Useful for detecting population-specific environmental effects
        GXE_ANALYSIS(
            ch_h2_results.collect { it[1] }
        )
        ch_gxe_results = GXE_ANALYSIS.out.gxe_estimates
        ch_versions = ch_versions.mix(GXE_ANALYSIS.out.versions)
    }

    // =========================================================================
    // CROSS-ANCESTRY GENETIC CORRELATION
    // =========================================================================

    if (genetic_correlation && h2_method == 'ldsc') {
        // Create pairs of ancestries for rg estimation
        ch_for_rg = LDSC_MUNGE.out.munged
            .map { meta, munged -> [[trait: meta.trait], meta.ancestry, munged] }
            .groupTuple(by: 0)
            .filter { trait_meta, ancestries, munged_files ->
                ancestries.size() >= 2
            }
            .flatMap { trait_meta, ancestries, munged_files ->
                def pairs = []
                for (int i = 0; i < ancestries.size() - 1; i++) {
                    for (int j = i + 1; j < ancestries.size(); j++) {
                        pairs << [
                            trait_meta + [anc1: ancestries[i], anc2: ancestries[j]],
                            munged_files[i],
                            munged_files[j]
                        ]
                    }
                }
                return pairs
            }

        LDSC_RG(
            ch_for_rg,
            ld_reference_dir
        )
        ch_rg_results = LDSC_RG.out.rg
        ch_versions = ch_versions.mix(LDSC_RG.out.versions)
    }

    // =========================================================================
    // COMPARE HERITABILITY ACROSS ANCESTRIES
    // =========================================================================

    H2_COMPARISON(
        ch_h2_results.collect { it[1] },
        ch_rg_results.collect { it[1] }.ifEmpty([]),
        ch_gxe_results.ifEmpty([])
    )
    ch_versions = ch_versions.mix(H2_COMPARISON.out.versions)

    emit:
    h2_estimates          = ch_h2_results                    // channel: [ meta, h2_results ]
    genetic_correlations  = ch_rg_results                    // channel: [ meta, rg_results ]
    gxe_estimates         = ch_gxe_results                   // channel: [ gxe_results ]
    comparison            = H2_COMPARISON.out.comparison     // channel: [ comparison_report ]
    versions              = ch_versions                      // channel: [ versions.yml ]
}
