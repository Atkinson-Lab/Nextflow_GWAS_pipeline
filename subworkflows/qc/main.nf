/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Quality control for genotype data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PLINK2_QC       } from '../../modules/local/plink2_qc'
include { PLINK2_CONVERT  } from '../../modules/local/plink2_convert'
include { QC_REPORT       } from '../../modules/local/qc_report'

workflow QC_WORKFLOW {
    take:
    ch_genotypes     // channel: [ meta, [ genotype_files ] ]
    maf_threshold    // value
    hwe_threshold    // value
    geno_missing     // value
    mind_missing     // value

    main:
    ch_versions = Channel.empty()

    // Convert to PLINK1 format if necessary for compatibility
    ch_to_convert = ch_genotypes.filter { meta, files ->
        meta.format != 'plink1'
    }
    ch_plink1 = ch_genotypes.filter { meta, files ->
        meta.format == 'plink1'
    }.map { meta, files ->
        [meta, files[0], files[1], files[2]]  // bed, bim, fam
    }

    if (ch_to_convert) {
        PLINK2_CONVERT(ch_to_convert)
        ch_plink1 = ch_plink1.mix(PLINK2_CONVERT.out.plink1)
        ch_versions = ch_versions.mix(PLINK2_CONVERT.out.versions)
    }

    // Run QC
    PLINK2_QC(
        ch_plink1,
        maf_threshold,
        hwe_threshold,
        geno_missing,
        mind_missing
    )
    ch_versions = ch_versions.mix(PLINK2_QC.out.versions)

    // Generate QC report
    QC_REPORT(
        PLINK2_QC.out.qc_stats.collect()
    )
    ch_versions = ch_versions.mix(QC_REPORT.out.versions)

    emit:
    genotypes = PLINK2_QC.out.genotypes  // channel: [ meta, bed, bim, fam ]
    qc_stats  = PLINK2_QC.out.qc_stats   // channel: [ meta, qc_files ]
    report    = QC_REPORT.out.report     // channel: [ report ]
    versions  = ch_versions              // channel: [ versions.yml ]
}
