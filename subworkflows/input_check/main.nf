/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_CHECK SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validates and parses input samplesheet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet  // file: path to samplesheet

    main:
    ch_versions = Channel.empty()

    SAMPLESHEET_CHECK(samplesheet)
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    // Parse validated samplesheet
    ch_genotypes = SAMPLESHEET_CHECK.out.csv
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id
            meta.cohort = row.cohort ?: 'default'

            def genotype_files = []
            if (row.bed && row.bim && row.fam) {
                // PLINK1 format
                genotype_files = [
                    file(row.bed, checkIfExists: true),
                    file(row.bim, checkIfExists: true),
                    file(row.fam, checkIfExists: true)
                ]
                meta.format = 'plink1'
            } else if (row.pgen && row.pvar && row.psam) {
                // PLINK2 format
                genotype_files = [
                    file(row.pgen, checkIfExists: true),
                    file(row.pvar, checkIfExists: true),
                    file(row.psam, checkIfExists: true)
                ]
                meta.format = 'plink2'
            } else if (row.bgen && row.sample) {
                // BGEN format
                genotype_files = [
                    file(row.bgen, checkIfExists: true),
                    file(row.sample, checkIfExists: true),
                    row.bgi ? file(row.bgi, checkIfExists: true) : null
                ].findAll { it != null }
                meta.format = 'bgen'
            } else if (row.vcf) {
                // VCF format
                genotype_files = [
                    file(row.vcf, checkIfExists: true),
                    row.vcf_index ? file(row.vcf_index, checkIfExists: true) : null
                ].findAll { it != null }
                meta.format = 'vcf'
            } else {
                error "Invalid samplesheet row: must provide genotype files"
            }

            [meta, genotype_files]
        }

    emit:
    genotypes = ch_genotypes  // channel: [ meta, [ files ] ]
    versions  = ch_versions   // channel: [ versions.yml ]
}
