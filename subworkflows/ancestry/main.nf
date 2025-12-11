/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANCESTRY_INFERENCE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry inference using GRAF-ANC and stratification
    Supports: EUR, AFR, EAS, SAS, AMR, MID, AAC, AHI, and HET (heterogeneous)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GRAF_ANCESTRY        } from '../../modules/local/graf_ancestry'
include { ADMIXTURE            } from '../../modules/local/admixture'
include { PCA_ANCESTRY         } from '../../modules/local/pca_ancestry'
include { ANCESTRY_STRATIFY    } from '../../modules/local/ancestry_stratify'
include { RFMIX_LOCAL_ANCESTRY } from '../../modules/local/rfmix'
include { GNOMIX_LOCAL_ANCESTRY } from '../../modules/local/gnomix'

workflow ANCESTRY_INFERENCE {
    take:
    ch_genotypes            // channel: [ meta, bed, bim, fam ]
    ancestry_method         // value: 'graf', 'admixture', 'pca'
    graf_reference          // file: GRAF reference
    run_local_ancestry      // boolean
    local_ancestry_method   // value: 'rfmix', 'gnomix', 'flare'

    main:
    ch_versions = Channel.empty()
    ch_ancestry_calls = Channel.empty()
    ch_local_ancestry = Channel.empty()

    // Run ancestry inference based on method
    if (ancestry_method == 'graf') {
        GRAF_ANCESTRY(
            ch_genotypes,
            graf_reference
        )
        ch_ancestry_calls = GRAF_ANCESTRY.out.ancestry_calls
        ch_versions = ch_versions.mix(GRAF_ANCESTRY.out.versions)
    } else if (ancestry_method == 'admixture') {
        ADMIXTURE(ch_genotypes)
        ch_ancestry_calls = ADMIXTURE.out.ancestry_calls
        ch_versions = ch_versions.mix(ADMIXTURE.out.versions)
    } else if (ancestry_method == 'pca') {
        PCA_ANCESTRY(ch_genotypes)
        ch_ancestry_calls = PCA_ANCESTRY.out.ancestry_calls
        ch_versions = ch_versions.mix(PCA_ANCESTRY.out.versions)
    }

    // Stratify genotypes by ancestry group
    ANCESTRY_STRATIFY(
        ch_genotypes.join(ch_ancestry_calls)
    )
    ch_stratified = ANCESTRY_STRATIFY.out.stratified_genotypes
    ch_versions = ch_versions.mix(ANCESTRY_STRATIFY.out.versions)

    // Run local ancestry inference if requested
    if (run_local_ancestry) {
        if (local_ancestry_method == 'rfmix') {
            RFMIX_LOCAL_ANCESTRY(
                ch_genotypes.join(ch_ancestry_calls)
            )
            ch_local_ancestry = RFMIX_LOCAL_ANCESTRY.out.local_ancestry
            ch_versions = ch_versions.mix(RFMIX_LOCAL_ANCESTRY.out.versions)
        } else if (local_ancestry_method == 'gnomix') {
            GNOMIX_LOCAL_ANCESTRY(
                ch_genotypes.join(ch_ancestry_calls)
            )
            ch_local_ancestry = GNOMIX_LOCAL_ANCESTRY.out.local_ancestry
            ch_versions = ch_versions.mix(GNOMIX_LOCAL_ANCESTRY.out.versions)
        }
    }

    emit:
    ancestry_calls       = ch_ancestry_calls         // channel: [ meta, ancestry_file ]
    stratified_genotypes = ch_stratified             // channel: [ meta, bed, bim, fam ] (meta includes ancestry)
    local_ancestry       = ch_local_ancestry         // channel: [ meta, local_ancestry_files ]
    versions             = ch_versions               // channel: [ versions.yml ]
}
