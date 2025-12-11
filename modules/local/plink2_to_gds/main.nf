process PLINK2_TO_GDS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-snprelate=1.36.0 bioconda::bioconductor-seqarray=1.42.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-snprelate:1.36.0' :
        'quay.io/biocontainers/bioconductor-snprelate:1.36.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${prefix}.gds"), emit: gds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry ?: 'all'}"
    """
    #!/usr/bin/env Rscript

    library(SNPRelate)
    library(SeqArray)

    # Convert PLINK to GDS format
    # GDS is required for GENESIS and provides efficient storage

    snpgdsBED2GDS(
        bed.fn = "${bed}",
        bim.fn = "${bim}",
        fam.fn = "${fam}",
        out.gdsfn = "${prefix}.gds",
        verbose = TRUE
    )

    # Verify GDS file
    gds <- snpgdsOpen("${prefix}.gds")
    cat("GDS file created successfully\\n")
    cat("Number of samples:", length(read.gdsn(index.gdsn(gds, "sample.id"))), "\\n")
    cat("Number of SNPs:", length(read.gdsn(index.gdsn(gds, "snp.id"))), "\\n")
    snpgdsClose(gds)

    # Write versions
    writeLines(c(
        '"${task.process}":',
        paste0('    snprelate: "', packageVersion("SNPRelate"), '"'),
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}
