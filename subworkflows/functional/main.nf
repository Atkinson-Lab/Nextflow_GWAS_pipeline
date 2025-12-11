/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONAL_ANNOT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Functional annotation and pathway analysis
    Tools: MAGMA, FUMA, LAVA, FLAMES
    Includes cell-type specific analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAGMA_ANNOTATE       } from '../../modules/local/magma_annotate'
include { MAGMA_GENE           } from '../../modules/local/magma_gene'
include { MAGMA_GENESET        } from '../../modules/local/magma_geneset'
include { MAGMA_CELLTYPE       } from '../../modules/local/magma_celltype'
include { FUMA_ANALYSIS        } from '../../modules/local/fuma_analysis'
include { LAVA                 } from '../../modules/local/lava'
include { FLAMES               } from '../../modules/local/flames'
include { CELL_TYPE_ENRICHMENT } from '../../modules/local/cell_type_enrichment'
include { FUNCTIONAL_SUMMARY   } from '../../modules/local/functional_summary'

workflow FUNCTIONAL_ANNOT {
    take:
    ch_gwas_results       // channel: [ meta, sumstats ]
    annotation_tools      // list: tools to run
    magma_geneset         // file: MAGMA gene set file
    magma_annotation      // file: MAGMA annotation file
    lava_partition        // file: LAVA partition file
    cell_type_analysis    // boolean
    cell_type_datasets    // file: cell-type expression datasets

    main:
    ch_versions = Channel.empty()
    ch_functional_results = Channel.empty()

    // =========================================================================
    // MAGMA ANALYSIS
    // =========================================================================

    if ('magma' in annotation_tools) {
        // Gene annotation
        MAGMA_ANNOTATE(
            ch_gwas_results,
            magma_annotation
        )
        ch_versions = ch_versions.mix(MAGMA_ANNOTATE.out.versions)

        // Gene-level analysis
        MAGMA_GENE(
            MAGMA_ANNOTATE.out.annotated
        )
        ch_functional_results = ch_functional_results.mix(
            MAGMA_GENE.out.genes.map { meta, results -> [meta + [tool: 'magma_gene'], results] }
        )
        ch_versions = ch_versions.mix(MAGMA_GENE.out.versions)

        // Gene-set enrichment analysis
        if (magma_geneset) {
            MAGMA_GENESET(
                MAGMA_GENE.out.genes,
                magma_geneset
            )
            ch_functional_results = ch_functional_results.mix(
                MAGMA_GENESET.out.geneset.map { meta, results -> [meta + [tool: 'magma_geneset'], results] }
            )
            ch_versions = ch_versions.mix(MAGMA_GENESET.out.versions)
        }

        // Cell-type specific analysis with MAGMA
        if (cell_type_analysis && cell_type_datasets) {
            MAGMA_CELLTYPE(
                MAGMA_GENE.out.genes,
                cell_type_datasets
            )
            ch_functional_results = ch_functional_results.mix(
                MAGMA_CELLTYPE.out.celltype.map { meta, results -> [meta + [tool: 'magma_celltype'], results] }
            )
            ch_versions = ch_versions.mix(MAGMA_CELLTYPE.out.versions)
        }
    }

    // =========================================================================
    // FUMA ANALYSIS
    // =========================================================================

    if ('fuma' in annotation_tools) {
        // Note: FUMA typically requires web interface, but we can prepare data
        // and potentially use FUMA API or local implementation
        FUMA_ANALYSIS(
            ch_gwas_results
        )
        ch_functional_results = ch_functional_results.mix(
            FUMA_ANALYSIS.out.results.map { meta, results -> [meta + [tool: 'fuma'], results] }
        )
        ch_versions = ch_versions.mix(FUMA_ANALYSIS.out.versions)
    }

    // =========================================================================
    // LAVA (Local Analysis of [co]Variant Association)
    // =========================================================================

    if ('lava' in annotation_tools) {
        LAVA(
            ch_gwas_results,
            lava_partition
        )
        ch_functional_results = ch_functional_results.mix(
            LAVA.out.results.map { meta, results -> [meta + [tool: 'lava'], results] }
        )
        ch_versions = ch_versions.mix(LAVA.out.versions)
    }

    // =========================================================================
    // FLAMES (Fine-mapping Linked Annotation of Molecular Effects at GWAS loci)
    // =========================================================================

    if ('flames' in annotation_tools) {
        FLAMES(
            ch_gwas_results
        )
        ch_functional_results = ch_functional_results.mix(
            FLAMES.out.results.map { meta, results -> [meta + [tool: 'flames'], results] }
        )
        ch_versions = ch_versions.mix(FLAMES.out.versions)
    }

    // =========================================================================
    // ADDITIONAL CELL-TYPE ENRICHMENT
    // =========================================================================

    if (cell_type_analysis && cell_type_datasets && !('magma' in annotation_tools)) {
        CELL_TYPE_ENRICHMENT(
            ch_gwas_results,
            cell_type_datasets
        )
        ch_functional_results = ch_functional_results.mix(
            CELL_TYPE_ENRICHMENT.out.results.map { meta, results -> [meta + [tool: 'celltype'], results] }
        )
        ch_versions = ch_versions.mix(CELL_TYPE_ENRICHMENT.out.versions)
    }

    // =========================================================================
    // SUMMARIZE FUNCTIONAL RESULTS
    // =========================================================================

    FUNCTIONAL_SUMMARY(
        ch_functional_results.collect { it[1] }
    )
    ch_versions = ch_versions.mix(FUNCTIONAL_SUMMARY.out.versions)

    emit:
    annotation_results = ch_functional_results              // channel: [ meta, results ]
    summary            = FUNCTIONAL_SUMMARY.out.summary     // channel: [ summary_file ]
    versions           = ch_versions                        // channel: [ versions.yml ]
}
