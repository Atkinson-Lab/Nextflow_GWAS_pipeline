/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORTING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generate comprehensive analysis reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC        } from '../../modules/local/multiqc'
include { CUSTOM_REPORT  } from '../../modules/local/custom_report'

workflow REPORTING {
    take:
    ch_multiqc_files    // channel: collected files for MultiQC
    ch_gwas_results     // channel: collected GWAS results
    ch_meta_results     // channel: collected meta-analysis results
    ch_fm_results       // channel: collected fine-mapping results
    ch_h2_results       // channel: collected heritability results
    ch_prs_comparison   // channel: collected PRS comparison results

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // MULTIQC REPORT
    // =========================================================================

    MULTIQC(
        ch_multiqc_files
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // =========================================================================
    // CUSTOM PIPELINE REPORT
    // =========================================================================

    // Generate comprehensive HTML report with all analysis summaries
    CUSTOM_REPORT(
        ch_gwas_results,
        ch_meta_results,
        ch_fm_results,
        ch_h2_results,
        ch_prs_comparison
    )
    ch_versions = ch_versions.mix(CUSTOM_REPORT.out.versions)

    emit:
    multiqc_report  = MULTIQC.out.report         // channel: [ multiqc_report ]
    custom_report   = CUSTOM_REPORT.out.report   // channel: [ custom_report ]
    versions        = ch_versions                // channel: [ versions.yml ]
}
