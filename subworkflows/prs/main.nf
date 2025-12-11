/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Multi-ancestry PRS calculation and validation
    Methods: PRS-CSx, PRS-CS, GAUDI, LDpred2, PRS-weighted
    Includes validation framework and best method determination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PRS_CSX           } from '../../modules/local/prs_csx'
include { PRS_CS            } from '../../modules/local/prs_cs'
include { GAUDI             } from '../../modules/local/gaudi'
include { LDPRED2           } from '../../modules/local/ldpred2'
include { PRS_WEIGHTED      } from '../../modules/local/prs_weighted'
include { PRS_CALCULATE     } from '../../modules/local/prs_calculate'
include { PRS_VALIDATION    } from '../../modules/local/prs_validation'
include { PRS_COMPARISON    } from '../../modules/local/prs_comparison'
include { PRS_BEST_METHOD   } from '../../modules/local/prs_best_method'

workflow PRS_WORKFLOW {
    take:
    ch_sumstats          // channel: [ meta, sumstats, source_type ]
    ch_genotypes         // channel: [ meta, bed, bim, fam ]
    ch_local_ancestry    // channel: [ meta, local_ancestry_files ]
    prs_methods          // list: methods to run
    prscsx_reference     // file: PRS-CSx LD reference
    ld_reference_dir     // file: LD reference directory
    validation_cohort    // file: validation cohort data
    run_validation       // boolean
    determine_best       // boolean

    main:
    ch_versions = Channel.empty()
    ch_prs_weights = Channel.empty()
    ch_prs_scores = Channel.empty()
    ch_validation = Channel.empty()

    // =========================================================================
    // CALCULATE PRS WEIGHTS USING MULTIPLE METHODS
    // =========================================================================

    // PRS-CSx (multi-ancestry, uses all populations jointly)
    if ('prs-csx' in prs_methods) {
        // Group sumstats by trait for multi-ancestry PRS-CSx
        ch_prscsx_input = ch_sumstats
            .filter { meta, sumstats, source -> source == 'meta' || source != 'meta' }
            .map { meta, sumstats, source -> [[trait: meta.trait], sumstats, meta.ancestry ?: 'meta'] }
            .groupTuple(by: 0)

        PRS_CSX(
            ch_prscsx_input,
            prscsx_reference
        )
        ch_prs_weights = ch_prs_weights.mix(
            PRS_CSX.out.weights.map { meta, weights -> [meta + [method: 'prs-csx'], weights] }
        )
        ch_versions = ch_versions.mix(PRS_CSX.out.versions)
    }

    // PRS-CS (single ancestry)
    if ('prs-cs' in prs_methods) {
        ch_prscs_input = ch_sumstats
            .filter { meta, sumstats, source -> source != 'meta' }
            .map { meta, sumstats, source -> [meta, sumstats] }

        PRS_CS(
            ch_prscs_input,
            ld_reference_dir
        )
        ch_prs_weights = ch_prs_weights.mix(
            PRS_CS.out.weights.map { meta, weights -> [meta + [method: 'prs-cs'], weights] }
        )
        ch_versions = ch_versions.mix(PRS_CS.out.versions)
    }

    // GAUDI (local ancestry-informed PRS)
    if ('gaudi' in prs_methods && ch_local_ancestry) {
        ch_gaudi_input = ch_sumstats
            .map { meta, sumstats, source -> [[trait: meta.trait], sumstats, source] }
            .groupTuple(by: 0)
            .combine(ch_local_ancestry)

        GAUDI(
            ch_gaudi_input,
            ld_reference_dir
        )
        ch_prs_weights = ch_prs_weights.mix(
            GAUDI.out.weights.map { meta, weights -> [meta + [method: 'gaudi'], weights] }
        )
        ch_versions = ch_versions.mix(GAUDI.out.versions)
    }

    // LDpred2 (single ancestry, Bayesian)
    if ('ldpred2' in prs_methods) {
        ch_ldpred2_input = ch_sumstats
            .filter { meta, sumstats, source -> source != 'meta' }
            .map { meta, sumstats, source -> [meta, sumstats] }

        LDPRED2(
            ch_ldpred2_input,
            ld_reference_dir
        )
        ch_prs_weights = ch_prs_weights.mix(
            LDPRED2.out.weights.map { meta, weights -> [meta + [method: 'ldpred2'], weights] }
        )
        ch_versions = ch_versions.mix(LDPRED2.out.versions)
    }

    // PRS-weighted (optimally weighted across ancestries)
    if ('prs-weighted' in prs_methods) {
        // Requires ancestry-specific PRS first
        ch_weights_for_weighted = ch_prs_weights
            .filter { meta, weights -> meta.method in ['prs-cs', 'ldpred2'] }
            .map { meta, weights -> [[trait: meta.trait], meta.ancestry, weights] }
            .groupTuple(by: 0)

        PRS_WEIGHTED(
            ch_weights_for_weighted
        )
        ch_prs_weights = ch_prs_weights.mix(
            PRS_WEIGHTED.out.weights.map { meta, weights -> [meta + [method: 'prs-weighted'], weights] }
        )
        ch_versions = ch_versions.mix(PRS_WEIGHTED.out.versions)
    }

    // =========================================================================
    // CALCULATE PRS SCORES
    // =========================================================================

    // Calculate PRS for each individual using each method's weights
    ch_calc_input = ch_prs_weights.combine(ch_genotypes)

    PRS_CALCULATE(
        ch_calc_input
    )
    ch_prs_scores = PRS_CALCULATE.out.scores
    ch_versions = ch_versions.mix(PRS_CALCULATE.out.versions)

    // =========================================================================
    // VALIDATION AND METHOD COMPARISON
    // =========================================================================

    if (run_validation && validation_cohort) {
        PRS_VALIDATION(
            ch_prs_scores,
            validation_cohort
        )
        ch_validation = PRS_VALIDATION.out.results
        ch_versions = ch_versions.mix(PRS_VALIDATION.out.versions)
    }

    // Compare methods and determine best approach
    if (determine_best) {
        // Group scores by trait for comparison
        ch_comparison_input = ch_prs_scores
            .map { meta, scores -> [[trait: meta.trait], meta.method, meta.ancestry ?: 'multi', scores] }
            .groupTuple(by: 0)

        PRS_COMPARISON(
            ch_comparison_input,
            ch_validation.collect { it[1] }.ifEmpty([])
        )
        ch_versions = ch_versions.mix(PRS_COMPARISON.out.versions)

        // Determine best method overall and per ancestry
        PRS_BEST_METHOD(
            PRS_COMPARISON.out.comparison
        )
        ch_best = PRS_BEST_METHOD.out.best_method
        ch_versions = ch_versions.mix(PRS_BEST_METHOD.out.versions)
    }

    emit:
    prs_weights        = ch_prs_weights                                           // channel: [ meta, weights ]
    prs_scores         = ch_prs_scores                                            // channel: [ meta, scores ]
    validation_results = ch_validation                                            // channel: [ meta, validation ]
    method_comparison  = determine_best ? PRS_COMPARISON.out.comparison : Channel.empty()  // channel: [ meta, comparison ]
    best_method        = determine_best ? ch_best : Channel.empty()               // channel: [ best_method_report ]
    versions           = ch_versions                                              // channel: [ versions.yml ]
}
