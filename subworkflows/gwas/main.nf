/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GWAS_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry-stratified GWAS using multiple tools
    Supports: REGENIE, SAIGE, BOLT-LMM, PLINK2, GENESIS, SPA-Cox
    Admixed-optimized: GENESIS, SAIGE, Tractor (local ancestry-aware)

    Tractor combines LAT1+LAT2 for Latino analysis (3-way: EUR,AFR,NAT)
    Tractor runs AAC separately (2-way: EUR,AFR)
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
include { PLINK_TO_VCF        } from '../../modules/local/format_conversion'
include { MERGE_PLINK_FILES   } from '../../modules/local/format_conversion'
include { VALIDATE_GWAS_INPUT } from '../../modules/local/validation'

workflow GWAS_WORKFLOW {
    take:
    ch_genotypes       // channel: [ meta, bed, bim, fam ] (meta includes ancestry, trait)
    ch_phenotypes      // channel: [ meta, phenotype_file ]
    covariate_cols     // value: covariate column names
    gwas_tool          // value: 'regenie', 'saige', 'bolt-lmm', 'plink2', 'genesis'
    gwas_model         // value: 'additive', 'dominant', 'recessive'
    kinship_matrix     // path: kinship/GRM matrix (optional, for GENESIS)
    run_tractor        // boolean: run Tractor for admixed populations
    ch_local_ancestry  // channel: [ meta, msp_file ] local ancestry calls (for Tractor)
    tractor_aac_pops   // value: ancestral populations for AAC (e.g., "EUR,AFR")
    tractor_lat_pops   // value: ancestral populations for Latino (e.g., "EUR,AFR,NAT")
    survival_analysis  // boolean: run survival/time-to-event analysis
    time_col           // value: time column for survival analysis
    event_col          // value: event column for survival analysis

    main:
    ch_versions = Channel.empty()
    ch_gwas_results = Channel.empty()
    ch_tractor_results = Channel.empty()
    ch_survival_results = Channel.empty()

    // =========================================================================
    // INPUT VALIDATION
    // =========================================================================

    VALIDATE_GWAS_INPUT(
        ch_genotypes,
        ch_phenotypes.map { meta, pheno -> pheno }.first()
    )
    ch_versions = ch_versions.mix(VALIDATE_GWAS_INPUT.out.versions)

    // Combine genotypes with phenotypes
    ch_gwas_input = ch_genotypes.combine(ch_phenotypes.map { meta, pheno -> pheno })

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
            .map { meta, bed, bim, fam, pheno ->
                [[id: meta.id, ancestry: meta.ancestry], meta, bed, bim, fam, pheno]
            }
            .combine(REGENIE_STEP1.out.predictions.map { meta, pred -> [meta, pred] }, by: 0)
            .map { key, meta, bed, bim, fam, pheno, predictions ->
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
                [[id: meta.id, ancestry: meta.ancestry, trait: meta.trait, binary: meta.binary ?: false],
                 bed, bim, fam, pheno]
            }

        SAIGE_STEP1(
            ch_step1_input,
            covariate_cols
        )
        ch_versions = ch_versions.mix(SAIGE_STEP1.out.versions)

        ch_step2_input = ch_step1_input
            .map { meta, bed, bim, fam, pheno -> [meta, bed, bim, fam] }
            .join(SAIGE_STEP1.out.model)

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

        // Prepare input for GENESIS null model
        ch_genesis_input = PLINK2_TO_GDS.out.gds
            .combine(ch_phenotypes.map { meta, pheno -> pheno })

        // Fit null model
        GENESIS_NULL_MODEL(
            ch_genesis_input,
            covariate_cols,
            kinship_matrix ?: []
        )
        ch_versions = ch_versions.mix(GENESIS_NULL_MODEL.out.versions)

        // Run association
        ch_assoc_input = PLINK2_TO_GDS.out.gds
            .join(GENESIS_NULL_MODEL.out.null_model)

        GENESIS_ASSOC(
            ch_assoc_input,
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
    // - AAC (African American): 2-way admixture (EUR, AFR)
    // - Latino (LAT1 + LAT2 COMBINED): 3-way admixture (EUR, AFR, NAT)
    //   LAT1 = Mexican/Central American, LAT2 = Caribbean/South American
    //   These are combined for better statistical power in Tractor analysis
    // Decomposes genetic effects by ancestral origin

    if (run_tractor && ch_local_ancestry) {
        // =====================================================================
        // AFRICAN AMERICAN (AAC) - 2-way admixture
        // =====================================================================
        ch_aac = ch_genotypes
            .filter { meta, bed, bim, fam ->
                meta.ancestry == 'AAC'
            }
            .join(ch_local_ancestry)

        // Convert PLINK to VCF for Tractor (AAC)
        PLINK_TO_VCF(
            ch_aac.map { meta, bed, bim, fam, la -> [meta, bed, bim, fam] },
            'aac'
        )

        ch_tractor_aac_input = PLINK_TO_VCF.out.vcf
            .join(ch_aac.map { meta, bed, bim, fam, la -> [meta, la] })
            .map { meta, vcf, vcf_idx, la_files ->
                [meta + [tractor_pops: tractor_aac_pops, tractor_group: 'AAC'], vcf, la_files]
            }

        // Extract ancestry-specific tract dosages (AAC)
        TRACTOR_EXTRACT_TRACTS(
            ch_tractor_aac_input,
            tractor_aac_pops
        )
        ch_versions = ch_versions.mix(TRACTOR_EXTRACT_TRACTS.out.versions.first())

        // Run Tractor association (AAC)
        ch_tractor_aac_assoc = TRACTOR_EXTRACT_TRACTS.out.ancestry_dosages
            .combine(ch_phenotypes.map { meta, pheno -> pheno })

        TRACTOR_ASSOC(
            ch_tractor_aac_assoc,
            covariate_cols
        )
        ch_tractor_results = ch_tractor_results.mix(TRACTOR_ASSOC.out.joint_results)
        ch_versions = ch_versions.mix(TRACTOR_ASSOC.out.versions.first())

        // =====================================================================
        // LATINO (LAT1 + LAT2 COMBINED) - 3-way admixture
        // =====================================================================
        // Combine LAT1 and LAT2 samples for Tractor analysis
        // This provides better statistical power for 3-way admixture modeling

        ch_latino_separate = ch_genotypes
            .filter { meta, bed, bim, fam ->
                meta.ancestry in ['LAT1', 'LAT2', 'AHI']
            }
            .join(ch_local_ancestry)

        // Merge LAT1 and LAT2 PLINK files into single "LATINO" group
        // Group by study ID (samples from same cohort get merged)
        ch_latino_for_merge = ch_latino_separate
            .map { meta, bed, bim, fam, la ->
                def merge_key = meta.id.replaceAll(/\.(LAT1|LAT2|AHI)$/, '')
                [[merge_id: merge_key, ancestry: 'LATINO'], meta, bed, bim, fam, la]
            }
            .groupTuple(by: 0)
            .map { merge_meta, metas, beds, bims, fams, las ->
                // Collect all files for merging
                [merge_meta, beds, bims, fams, las, metas.collect { it.ancestry }]
            }

        // Check if we have data to merge
        ch_latino_for_merge
            .branch {
                single: it[1].size() == 1
                multiple: it[1].size() > 1
            }
            .set { ch_latino_branched }

        // For single ancestry (no merge needed)
        ch_latino_single = ch_latino_branched.single
            .map { merge_meta, beds, bims, fams, las, ancestries ->
                def new_meta = [
                    id: "${merge_meta.merge_id}.LATINO",
                    ancestry: 'LATINO',
                    original_ancestries: ancestries.join(','),
                    n_samples: 0  // Will be computed
                ]
                [new_meta, beds[0], bims[0], fams[0], las[0]]
            }

        // For multiple ancestries - merge PLINK files
        MERGE_PLINK_FILES(
            ch_latino_branched.multiple
                .map { merge_meta, beds, bims, fams, las, ancestries ->
                    [[id: "${merge_meta.merge_id}.LATINO", ancestry: 'LATINO',
                      original_ancestries: ancestries.join(',')],
                     beds, bims, fams]
                }
        )

        // Combine merged local ancestry files for Latino
        ch_latino_merged_la = ch_latino_branched.multiple
            .map { merge_meta, beds, bims, fams, las, ancestries ->
                [[id: "${merge_meta.merge_id}.LATINO"], las]
            }

        // Combine single and merged Latino data
        ch_latino_combined = ch_latino_single
            .mix(
                MERGE_PLINK_FILES.out.merged
                    .join(ch_latino_merged_la)
                    .map { meta, bed, bim, fam, las ->
                        // Use first LA file (should be combined upstream if needed)
                        [meta, bed, bim, fam, las[0]]
                    }
            )

        // Convert PLINK to VCF for Tractor (Latino)
        PLINK_TO_VCF(
            ch_latino_combined.map { meta, bed, bim, fam, la -> [meta, bed, bim, fam] },
            'latino'
        )

        ch_tractor_latino_input = PLINK_TO_VCF.out.vcf
            .join(ch_latino_combined.map { meta, bed, bim, fam, la -> [meta, la] })
            .map { meta, vcf, vcf_idx, la_files ->
                [meta + [tractor_pops: tractor_lat_pops, tractor_group: 'LATINO'], vcf, la_files]
            }

        // Extract ancestry-specific tract dosages (Latino)
        TRACTOR_EXTRACT_TRACTS(
            ch_tractor_latino_input,
            tractor_lat_pops
        )

        // Run Tractor association (Latino)
        ch_tractor_latino_assoc = TRACTOR_EXTRACT_TRACTS.out.ancestry_dosages
            .combine(ch_phenotypes.map { meta, pheno -> pheno })

        TRACTOR_ASSOC(
            ch_tractor_latino_assoc,
            covariate_cols
        )
        ch_tractor_results = ch_tractor_results.mix(TRACTOR_ASSOC.out.joint_results)
    }

    // =========================================================================
    // SURVIVAL / TIME-TO-EVENT ANALYSIS (SPA-Cox)
    // =========================================================================

    if (survival_analysis && time_col && event_col) {
        // Filter to traits with time-to-event data
        ch_survival_input = ch_gwas_input

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

    // Only run filtering if we have results
    ch_all_gwas
        .ifEmpty { log.warn "No GWAS results to filter" }
        .set { ch_gwas_to_filter }

    GWAS_FILTER(
        ch_gwas_to_filter
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
