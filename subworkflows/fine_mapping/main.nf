/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINE_MAPPING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Within-ancestry and multi-ancestry fine-mapping
    Within-ancestry: PolyFun+SuSIE, SuSIE, FINEMAP, CAVIAR
    Multi-ancestry: MG-FLASH-FM (related traits), SuSIE-ME (single/unrelated traits)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { POLYFUN_MUNGING    } from '../../modules/local/polyfun_munging'
include { POLYFUN_SUSIE      } from '../../modules/local/polyfun_susie'
include { SUSIE              } from '../../modules/local/susie'
include { FINEMAP            } from '../../modules/local/finemap'
include { CAVIAR             } from '../../modules/local/caviar'
include { MG_FLASH_FM        } from '../../modules/local/mg_flash_fm'
include { SUSIE_ME           } from '../../modules/local/susie_me'
include { CREDIBLE_SETS      } from '../../modules/local/credible_sets'
include { EXTRACT_REGIONS    } from '../../modules/local/extract_regions'

workflow FINE_MAPPING {
    take:
    ch_within_ancestry_input   // channel: [ meta, sumstats, method ]
    ch_multi_ancestry_input    // channel: [ [trait], [ancestries], [sumstats_files] ]
    multi_ancestry_method      // value: 'mg-flash-fm' or 'susie-me'
    ld_reference_dir           // file: directory with ancestry-specific LD panels
    polyfun_annotations        // file: PolyFun annotation files
    window_kb                  // value: fine-mapping window in kb
    max_causal                 // value: maximum causal variants

    main:
    ch_versions = Channel.empty()
    ch_within_results = Channel.empty()
    ch_multi_results = Channel.empty()
    ch_credible_sets = Channel.empty()

    // =========================================================================
    // WITHIN-ANCESTRY FINE-MAPPING
    // =========================================================================

    // Extract significant regions for fine-mapping
    ch_to_finemmap = ch_within_ancestry_input
        .map { meta, sumstats, method -> [meta, sumstats] }

    EXTRACT_REGIONS(
        ch_to_finemmap,
        window_kb
    )
    ch_versions = ch_versions.mix(EXTRACT_REGIONS.out.versions)

    // Join with method info
    ch_regions_with_method = EXTRACT_REGIONS.out.regions
        .join(ch_within_ancestry_input.map { meta, sumstats, method -> [meta, method] })

    // Branch by fine-mapping method
    ch_regions_with_method.branch {
        polyfun_susie: it[2] == 'polyfun-susie'
        susie: it[2] == 'susie'
        finemap: it[2] == 'finemap'
        caviar: it[2] == 'caviar'
    }.set { ch_fm_branched }

    // PolyFun+SuSIE (annotation-prioritized)
    if (polyfun_annotations) {
        // First munge sumstats for PolyFun
        POLYFUN_MUNGING(
            ch_fm_branched.polyfun_susie.map { meta, regions, method -> [meta, regions] },
            polyfun_annotations
        )
        ch_versions = ch_versions.mix(POLYFUN_MUNGING.out.versions)

        POLYFUN_SUSIE(
            POLYFUN_MUNGING.out.munged,
            ld_reference_dir,
            max_causal
        )
        ch_within_results = ch_within_results.mix(POLYFUN_SUSIE.out.results)
        ch_versions = ch_versions.mix(POLYFUN_SUSIE.out.versions)
    }

    // Standard SuSIE
    SUSIE(
        ch_fm_branched.susie.map { meta, regions, method -> [meta, regions] },
        ld_reference_dir,
        max_causal
    )
    ch_within_results = ch_within_results.mix(SUSIE.out.results)
    ch_versions = ch_versions.mix(SUSIE.out.versions)

    // FINEMAP
    FINEMAP(
        ch_fm_branched.finemap.map { meta, regions, method -> [meta, regions] },
        ld_reference_dir,
        max_causal
    )
    ch_within_results = ch_within_results.mix(FINEMAP.out.results)
    ch_versions = ch_versions.mix(FINEMAP.out.versions)

    // CAVIAR
    CAVIAR(
        ch_fm_branched.caviar.map { meta, regions, method -> [meta, regions] },
        ld_reference_dir,
        max_causal
    )
    ch_within_results = ch_within_results.mix(CAVIAR.out.results)
    ch_versions = ch_versions.mix(CAVIAR.out.versions)

    // =========================================================================
    // MULTI-ANCESTRY FINE-MAPPING
    // =========================================================================

    // Filter to traits with fine-mapping results from multiple ancestries
    ch_multi_fm_input = ch_multi_ancestry_input
        .filter { trait_meta, ancestries, sumstats ->
            ancestries.size() >= 2
        }

    if (multi_ancestry_method == 'mg-flash-fm') {
        // MG-FLASH-FM for related traits (multi-trait multi-ancestry)
        MG_FLASH_FM(
            ch_multi_fm_input,
            ld_reference_dir,
            window_kb,
            max_causal
        )
        ch_multi_results = MG_FLASH_FM.out.results
        ch_versions = ch_versions.mix(MG_FLASH_FM.out.versions)

    } else if (multi_ancestry_method == 'susie-me') {
        // SuSIE-ME for single/unrelated traits
        SUSIE_ME(
            ch_multi_fm_input,
            ld_reference_dir,
            max_causal
        )
        ch_multi_results = SUSIE_ME.out.results
        ch_versions = ch_versions.mix(SUSIE_ME.out.versions)
    }

    // =========================================================================
    // GENERATE CREDIBLE SETS
    // =========================================================================

    // Combine within and multi-ancestry results
    ch_all_fm_results = ch_within_results.mix(ch_multi_results)

    CREDIBLE_SETS(
        ch_all_fm_results
    )
    ch_credible_sets = CREDIBLE_SETS.out.credible_sets
    ch_versions = ch_versions.mix(CREDIBLE_SETS.out.versions)

    emit:
    fine_mapping_results = ch_all_fm_results        // channel: [ meta, fm_results ]
    credible_sets        = ch_credible_sets         // channel: [ meta, credible_set_file ]
    within_ancestry      = ch_within_results        // channel: [ meta, within_results ]
    multi_ancestry       = ch_multi_results         // channel: [ meta, multi_results ]
    versions             = ch_versions              // channel: [ versions.yml ]
}
