process BOLT_LMM {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    conda "bioconda::bolt-lmm=2.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bolt-lmm:2.4.1--h9f5acd7_1' :
        'quay.io/biocontainers/bolt-lmm:2.4.1--h9f5acd7_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype)
    val covariate_cols
    val gwas_model

    output:
    tuple val(meta), path("${prefix}.bolt.stats.gz"), emit: summary_stats
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ? "--covarCol=${covariate_cols.split(',').join(' --covarCol=')}" : ""
    """
    bolt \\
        --bfile=${bed.baseName} \\
        --phenoFile=${phenotype} \\
        --phenoCol=${meta.trait} \\
        ${covar_arg} \\
        --lmm \\
        --LDscoresUseChip \\
        --statsFile=${prefix}.bolt.stats \\
        --numThreads=${task.cpus} \\
        $args \\
        2>&1 | tee ${prefix}.log

    gzip ${prefix}.bolt.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt-lmm: \$(bolt --version 2>&1 | head -1 | sed 's/BOLT-LMM v//')
    END_VERSIONS
    """
}
