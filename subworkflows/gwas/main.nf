/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GWAS_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry-stratified GWAS using multiple tools
    Supports: REGENIE, SAIGE, BOLT-LMM, PLINK2, GENESIS, SPA-Cox
    Admixed-optimized: GENESIS, SAIGE, Tractor (local ancestry-aware)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REGENIE_STEP1       } from '../../modules/local/regenie_step1'
include { REGENIE_STEP2       } from '../../modules/local/regenie_step2'
include { SAIGE_STEP1         } from '../../modules/local/saige_step1'
include { SAIGE_STEP2         } from '../../modules/local/saige_step2'
include { BOLT_LMM            } from '../../modules/local/bolt_lmm'
include { PLINK2_GWAS         } from '../../modules/local/plink2_gwas'
include { GENESIS_NULL_MODEL  } from '../../modules/local/genesis'
include { GENESIS_ASSOC       } from '../../modules/local/genesis'
include { SPA_COX_STEP1       } from '../../modules/local/spa_cox'
include { SPA_COX_STEP2       } from '../../modules/local/spa_cox'
include { TRACTOR_EXTRACT_TRACTS } from '../../modules/local/tractor'
include { TRACTOR_ASSOC       } from '../../modules/local/tractor'
include { GWAS_FILTER         } from '../../modules/local/gwas_filter'
include { CLUMP_REGIONS       } from '../../modules/local/clump_regions'
include { PLINK2_TO_GDS       } from '../../modules/local/plink2_to_gds'

workflow GWAS_WORKFLOW {
    take:
    ch_genotypes       // channel: [ meta, bed, bim, fam ] (meta includes ancestry, trait)
    ch_phenotypes      // channel: [ phenotype_file ]
    covariate_cols     // value: covariate column names
    gwas_tool          // value: 'regenie', 'saige', 'bolt-lmm', 'plink2', 'genesis'
    gwas_model         // value: 'additive', 'dominant', 'recessive'
    kinship_matrix     // path: kinship/GRM matrix (optional, for GENESIS)
    run_tractor        // boolean: run Tractor for admixed populations
    ch_local_ancestry  // channel: local ancestry calls (for Tractor)
    tractor_ancestries // value: ancestral populations for Tractor (e.g., "EUR,AFR" or "EUR,AFR,NAT")
    survival_analysis  // boolean: run survival/time-to-event analysis
    time_col           // value: time column for survival analysis
    event_col          // value: event column for survival analysis

    main:
    ch_versions = Channel.empty()
    ch_gwas_results = Channel.empty()
    ch_tractor_results = Channel.empty()
    ch_survival_results = Channel.empty()

    // Combine genotypes with phenotypes
    ch_gwas_input = ch_genotypes.combine(ch_phenotypes)

    // =========================================================================
    // STANDARD GWAS TOOLS
    // =========================================================================

    if (gwas_tool == 'regenie') {
        // REGENIE Step 1: Fit null model
        // Group by ancestry for Step 1 (fit once per ancestry)
        ch_step1_input = ch_gwas_input
            .map { meta, bed, bim, fam, pheno ->
                [[id: meta.id, ancestry: meta.ancestry], bed, bim, fam, pheno]
            }
            .groupTuple(by: 0)
            .map { meta, beds, bims, fams, phenos ->
                [meta, beds[0], bims[0], fams[0], phenos[0]]
            }

        REGENIE_STEP1(
            ch_step1_input,
            covariate_cols
        )
        ch_versions = ch_versions.mix(REGENIE_STEP1.out.versions)

        // REGENIE Step 2: Association testing (per trait)
        ch_step2_input = ch_gwas_input
            .combine(REGENIE_STEP1.out.predictions, by: { a, b ->
                a[0].ancestry == b[0].ancestry
            })
            .map { meta, bed, bim, fam, pheno, pred_meta, predictions ->
                [meta, bed, bim, fam, pheno, predictions]
            }

        REGENIE_STEP2(
            ch_step2_input,
            covariate_cols,
            gwas_model
        )
        ch_gwas_results = REGENIE_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(REGENIE_STEP2.out.versions)

    } else if (gwas_tool == 'saige') {
        // SAIGE - optimized for admixed populations and case-control imbalance
        ch_step1_input = ch_gwas_input
            .map { meta, bed, bim, fam, pheno ->
                [[id: meta.id, ancestry: meta.ancestry, trait: meta.trait, binary: meta.binary],
                 bed, bim, fam, pheno]
            }

        SAIGE_STEP1(
            ch_step1_input,
            covariate_cols
        )
        ch_versions = ch_versions.mix(SAIGE_STEP1.out.versions)

        ch_step2_input = ch_gwas_input
            .join(SAIGE_STEP1.out.model, by: 0)

        SAIGE_STEP2(
            ch_step2_input,
            gwas_model
        )
        ch_gwas_results = SAIGE_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(SAIGE_STEP2.out.versions)

    } else if (gwas_tool == 'genesis') {
        // GENESIS - optimized for admixed populations with complex relatedness
        // Requires GDS format and optionally kinship matrix

        // Convert PLINK to GDS
        PLINK2_TO_GDS(
            ch_genotypes
        )
        ch_versions = ch_versions.mix(PLINK2_TO_GDS.out.versions)

        // Fit null model
        GENESIS_NULL_MODEL(
            PLINK2_TO_GDS.out.gds.join(ch_phenotypes),
            covariate_cols,
            kinship_matrix ?: []
        )
        ch_versions = ch_versions.mix(GENESIS_NULL_MODEL.out.versions)

        // Run association
        GENESIS_ASSOC(
            PLINK2_TO_GDS.out.gds.join(GENESIS_NULL_MODEL.out.null_model),
            'Score'  // Test type: Score, Wald, BinomiRare
        )
        ch_gwas_results = GENESIS_ASSOC.out.summary_stats
        ch_versions = ch_versions.mix(GENESIS_ASSOC.out.versions)

    } else if (gwas_tool == 'bolt-lmm') {
        BOLT_LMM(
            ch_gwas_input,
            covariate_cols,
            gwas_model
        )
        ch_gwas_results = BOLT_LMM.out.summary_stats
        ch_versions = ch_versions.mix(BOLT_LMM.out.versions)

    } else if (gwas_tool == 'plink2') {
        PLINK2_GWAS(
            ch_gwas_input,
            covariate_cols,
            gwas_model
        )
        ch_gwas_results = PLINK2_GWAS.out.summary_stats
        ch_versions = ch_versions.mix(PLINK2_GWAS.out.versions)
    }

    // =========================================================================
    // TRACTOR - Local Ancestry-Aware GWAS for Admixed Populations
    // =========================================================================
    // Run for: AAC (African American), LAT1 (Latino Type 1), LAT2 (Latino Type 2)
    // Decomposes genetic effects by ancestral origin

    if (run_tractor && ch_local_ancestry) {
        // Filter to admixed populations
        ch_admixed = ch_genotypes
            .filter { meta, bed, bim, fam ->
                meta.ancestry in ['AAC', 'AHI', 'LAT1', 'LAT2', 'ADMIXED']
            }
            .join(ch_local_ancestry)

        // Determine ancestral populations based on admixed group
        ch_tractor_input = ch_admixed
            .map { meta, bed, bim, fam, la_files ->
                def anc_pops = tractor_ancestries
                if (!anc_pops) {
                    // Default ancestral populations by admixed group
                    if (meta.ancestry == 'AAC') {
                        anc_pops = 'EUR,AFR'
                    } else if (meta.ancestry in ['AHI', 'LAT1', 'LAT2']) {
                        anc_pops = 'EUR,AFR,NAT'
                    } else {
                        anc_pops = 'EUR,AFR'  // Default 2-way
                    }
                }
                [meta + [tractor_pops: anc_pops], bed, bim, fam, la_files]
            }

        // Extract ancestry-specific tract dosages
        TRACTOR_EXTRACT_TRACTS(
            ch_tractor_input.map { meta, bed, bim, fam, la -> [meta, bed, la] },
            ch_tractor_input.map { it[0].tractor_pops }.first()
        )
        ch_versions = ch_versions.mix(TRACTOR_EXTRACT_TRACTS.out.versions)

        // Run Tractor association
        TRACTOR_ASSOC(
            TRACTOR_EXTRACT_TRACTS.out.ancestry_dosages.join(ch_phenotypes),
            covariate_cols
        )
        ch_tractor_results = TRACTOR_ASSOC.out.joint_results
        ch_versions = ch_versions.mix(TRACTOR_ASSOC.out.versions)
    }

    // =========================================================================
    // SURVIVAL / TIME-TO-EVENT ANALYSIS (SPA-Cox)
    // =========================================================================

    if (survival_analysis && time_col && event_col) {
        // Filter to traits with time-to-event data
        ch_survival_input = ch_gwas_input
            .filter { meta, bed, bim, fam, pheno ->
                // Check if phenotype file has required columns
                true  // Actual check would parse pheno file
            }

        SPA_COX_STEP1(
            ch_survival_input,
            covariate_cols,
            time_col,
            event_col
        )
        ch_versions = ch_versions.mix(SPA_COX_STEP1.out.versions)

        ch_cox_step2_input = ch_survival_input
            .map { meta, bed, bim, fam, pheno -> [meta, bed, bim, fam] }
            .join(SPA_COX_STEP1.out.model)

        SPA_COX_STEP2(
            ch_cox_step2_input
        )
        ch_survival_results = SPA_COX_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(SPA_COX_STEP2.out.versions)
    }

    // =========================================================================
    // FILTER AND CLUMP RESULTS
    // =========================================================================

    // Combine all GWAS results
    ch_all_gwas = ch_gwas_results
        .mix(ch_tractor_results)
        .mix(ch_survival_results)

    GWAS_FILTER(
        ch_all_gwas
    )
    ch_versions = ch_versions.mix(GWAS_FILTER.out.versions)

    // Identify independent signals by clumping
    CLUMP_REGIONS(
        GWAS_FILTER.out.significant,
        ch_genotypes.map { meta, bed, bim, fam -> [meta.ancestry, bed, bim, fam] }.unique()
    )
    ch_versions = ch_versions.mix(CLUMP_REGIONS.out.versions)

    emit:
    summary_stats     = ch_gwas_results                    // channel: [ meta, sumstats ]
    tractor_results   = ch_tractor_results                 // channel: [ meta, tractor_sumstats ]
    survival_results  = ch_survival_results                // channel: [ meta, survival_sumstats ]
    filtered_results  = GWAS_FILTER.out.filtered           // channel: [ meta, filtered_sumstats ]
    significant       = GWAS_FILTER.out.significant        // channel: [ meta, sig_variants ]
    clumped           = CLUMP_REGIONS.out.clumped          // channel: [ meta, clumped_results ]
    versions          = ch_versions                        // channel: [ versions.yml ]
}
