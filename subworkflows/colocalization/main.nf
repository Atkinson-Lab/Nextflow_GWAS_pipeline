/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLOCALIZATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Colocalization analysis with eQTL/sQTL/pQTL datasets
    Supports pooling of diverse QTL datasets for improved representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_QTL_DATA   } from '../../modules/local/prepare_qtl_data'
include { POOL_QTL_DATASETS  } from '../../modules/local/pool_qtl_datasets'
include { COLOC              } from '../../modules/local/coloc'
include { ECAVIAR            } from '../../modules/local/ecaviar'
include { FASTENLOC          } from '../../modules/local/fastenloc'
include { COLOC_SUMMARY      } from '../../modules/local/coloc_summary'

workflow COLOCALIZATION {
    take:
    ch_gwas_results    // channel: [ meta, sumstats ]
    ch_qtl_data        // channel: [ name, type, qtl_file ]
    coloc_method       // value: 'coloc', 'ecaviar', 'fastenloc'
    p1                 // value: prior for trait association
    p2                 // value: prior for QTL association
    p12                // value: prior for colocalization

    main:
    ch_versions = Channel.empty()
    ch_coloc_results = Channel.empty()

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
    // RUN COLOCALIZATION
    // =========================================================================

    // Create GWAS x QTL combinations for colocalization
    ch_coloc_input = ch_gwas_results.combine(ch_qtl_pooled)

    if (coloc_method == 'coloc') {
        COLOC(
            ch_coloc_input,
            p1,
            p2,
            p12
        )
        ch_coloc_results = COLOC.out.results
        ch_versions = ch_versions.mix(COLOC.out.versions)

    } else if (coloc_method == 'ecaviar') {
        ECAVIAR(
            ch_coloc_input
        )
        ch_coloc_results = ECAVIAR.out.results
        ch_versions = ch_versions.mix(ECAVIAR.out.versions)

    } else if (coloc_method == 'fastenloc') {
        FASTENLOC(
            ch_coloc_input
        )
        ch_coloc_results = FASTENLOC.out.results
        ch_versions = ch_versions.mix(FASTENLOC.out.versions)
    }

    // =========================================================================
    // SUMMARIZE COLOCALIZATION RESULTS
    // =========================================================================

    COLOC_SUMMARY(
        ch_coloc_results.collect { it[1] }
    )
    ch_versions = ch_versions.mix(COLOC_SUMMARY.out.versions)

    emit:
    coloc_results   = ch_coloc_results           // channel: [ meta, coloc_results ]
    summary         = COLOC_SUMMARY.out.summary  // channel: [ summary_file ]
    pooled_qtl      = ch_qtl_pooled              // channel: [ type, pooled_qtl_file ]
    versions        = ch_versions                // channel: [ versions.yml ]
}
