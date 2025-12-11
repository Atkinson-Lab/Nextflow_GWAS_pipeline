process LDSC_MUNGE {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_low'

    conda "bioconda::ldsc=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--pyhdfd78af_2' :
        'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(sumstats)

    output:
    tuple val(meta), path("${prefix}.sumstats.gz"), emit: munged
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    # Munge summary statistics for LDSC
    munge_sumstats.py \\
        --sumstats ${sumstats} \\
        --out ${prefix} \\
        --merge-alleles /opt/ldsc/w_hm3.snplist \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ldsc: \$(ldsc.py --version 2>&1 | head -1 || echo "1.0.1")
    END_VERSIONS
    """
}
