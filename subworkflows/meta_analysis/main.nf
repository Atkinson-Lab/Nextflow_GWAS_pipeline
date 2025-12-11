/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    META_ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Multi-ancestry meta-analysis using MR-MEGA
    Supports: MR-MEGA, MR-MEGA-env, METAL, METASOft
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MR_MEGA           } from '../../modules/local/mr_mega'
include { MR_MEGA_ENV       } from '../../modules/local/mr_mega_env'
include { METAL             } from '../../modules/local/metal'
include { METASOFT          } from '../../modules/local/metasoft'
include { HETEROGENEITY     } from '../../modules/local/heterogeneity'
include { MERGE_META_RESULTS } from '../../modules/local/merge_meta_results'

workflow META_ANALYSIS {
    take:
    ch_gwas_grouped    // channel: [ [trait, binary], [ancestries], [sumstats_files] ]
    meta_method        // value: 'mr-mega', 'metal', 'metasoft'
    run_mr_mega_env    // boolean
    env_variable       // value: environmental variable column

    main:
    ch_versions = Channel.empty()
    ch_meta_results = Channel.empty()
    ch_heterogeneity = Channel.empty()

    // Filter to traits with multiple ancestry groups
    ch_multi_ancestry = ch_gwas_grouped
        .filter { trait_meta, ancestries, sumstats ->
            ancestries.size() >= 2
        }

    if (meta_method == 'mr-mega') {
        // Standard MR-MEGA
        MR_MEGA(
            ch_multi_ancestry
        )
        ch_meta_results = MR_MEGA.out.meta_results
        ch_heterogeneity = MR_MEGA.out.heterogeneity
        ch_versions = ch_versions.mix(MR_MEGA.out.versions)

        // Optional: MR-MEGA with environmental variable
        if (run_mr_mega_env && env_variable) {
            MR_MEGA_ENV(
                ch_multi_ancestry,
                env_variable
            )
            // Add env results to output
            ch_meta_env_results = MR_MEGA_ENV.out.meta_results
            ch_versions = ch_versions.mix(MR_MEGA_ENV.out.versions)
        }

    } else if (meta_method == 'metal') {
        METAL(
            ch_multi_ancestry
        )
        ch_meta_results = METAL.out.meta_results
        ch_versions = ch_versions.mix(METAL.out.versions)

    } else if (meta_method == 'metasoft') {
        METASOFT(
            ch_multi_ancestry
        )
        ch_meta_results = METASOFT.out.meta_results
        ch_heterogeneity = METASOFT.out.heterogeneity
        ch_versions = ch_versions.mix(METASOFT.out.versions)
    }

    // Calculate heterogeneity statistics (I², Q)
    HETEROGENEITY(
        ch_meta_results
    )
    ch_versions = ch_versions.mix(HETEROGENEITY.out.versions)

    // Merge and format final meta-analysis results
    MERGE_META_RESULTS(
        ch_meta_results.join(HETEROGENEITY.out.stats)
    )
    ch_versions = ch_versions.mix(MERGE_META_RESULTS.out.versions)

    emit:
    meta_results     = MERGE_META_RESULTS.out.merged  // channel: [ meta, meta_sumstats ]
    heterogeneity    = HETEROGENEITY.out.stats        // channel: [ meta, het_stats ]
    ancestry_effects = ch_heterogeneity               // channel: [ meta, ancestry_effects ]
    versions         = ch_versions                    // channel: [ versions.yml ]
}
