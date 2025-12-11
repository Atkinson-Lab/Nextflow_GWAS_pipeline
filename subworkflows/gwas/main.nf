/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GWAS_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry-stratified GWAS using multiple tools
    Supports: REGENIE, SAIGE, BOLT-LMM, PLINK2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REGENIE_STEP1    } from '../../modules/local/regenie_step1'
include { REGENIE_STEP2    } from '../../modules/local/regenie_step2'
include { SAIGE_STEP1      } from '../../modules/local/saige_step1'
include { SAIGE_STEP2      } from '../../modules/local/saige_step2'
include { BOLT_LMM         } from '../../modules/local/bolt_lmm'
include { PLINK2_GWAS      } from '../../modules/local/plink2_gwas'
include { GWAS_FILTER      } from '../../modules/local/gwas_filter'
include { CLUMP_REGIONS    } from '../../modules/local/clump_regions'

workflow GWAS_WORKFLOW {
    take:
    ch_genotypes      // channel: [ meta, bed, bim, fam ] (meta includes ancestry, trait)
    ch_phenotypes     // channel: [ phenotype_file ]
    covariate_cols    // value: covariate column names
    gwas_tool         // value: 'regenie', 'saige', 'bolt-lmm', 'plink2'
    gwas_model        // value: 'additive', 'dominant', 'recessive'

    main:
    ch_versions = Channel.empty()
    ch_gwas_results = Channel.empty()

    // Combine genotypes with phenotypes
    ch_gwas_input = ch_genotypes.combine(ch_phenotypes)

    // Run GWAS based on tool selection
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
        // SAIGE Step 1: Fit null model
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

        // SAIGE Step 2: Association testing
        ch_step2_input = ch_gwas_input
            .join(SAIGE_STEP1.out.model, by: 0)

        SAIGE_STEP2(
            ch_step2_input,
            gwas_model
        )
        ch_gwas_results = SAIGE_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(SAIGE_STEP2.out.versions)

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

    // Filter and clump significant results
    GWAS_FILTER(
        ch_gwas_results
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
    filtered_results  = GWAS_FILTER.out.filtered           // channel: [ meta, filtered_sumstats ]
    significant       = GWAS_FILTER.out.significant        // channel: [ meta, sig_variants ]
    clumped           = CLUMP_REGIONS.out.clumped          // channel: [ meta, clumped_results ]
    versions          = ch_versions                        // channel: [ versions.yml ]
}
