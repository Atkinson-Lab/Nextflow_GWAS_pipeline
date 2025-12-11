process MAGMA_GENE {
    tag "$meta.trait"
    label 'process_medium'

    conda "bioconda::magma=1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magma:1.10--h9f5acd7_0' :
        'quay.io/biocontainers/magma:1.10--h9f5acd7_0' }"

    input:
    tuple val(meta), path(annotated_sumstats)

    output:
    tuple val(meta), path("${prefix}.genes.out"), emit: genes
    tuple val(meta), path("${prefix}.genes.raw"), emit: raw
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}"
    """
    # MAGMA gene-level analysis
    # Aggregates SNP-level associations to gene-level

    magma \\
        --bfile /opt/magma/g1000_eur \\
        --pval ${annotated_sumstats} N=\${N_SAMPLES:-10000} \\
        --gene-annot ${annotated_sumstats.baseName}.genes.annot \\
        --out ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma --version 2>&1 | head -1 | sed 's/MAGMA version: //')
    END_VERSIONS
    """
}
