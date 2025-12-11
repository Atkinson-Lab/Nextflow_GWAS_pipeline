process CLUMP_REGIONS {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_low'

    conda "bioconda::plink=1.90b6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h031d066_5' :
        'quay.io/biocontainers/plink:1.90b6.21--h031d066_5' }"

    input:
    tuple val(meta), path(significant)
    tuple val(ancestry), path(ld_bed), path(ld_bim), path(ld_fam)

    output:
    tuple val(meta), path("${prefix}.clumped"), emit: clumped
    tuple val(meta), path("${prefix}.clumped.ranges"), emit: ranges
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    # Clump significant variants to identify independent signals
    plink \\
        --bfile ${ld_bed.baseName} \\
        --clump ${significant} \\
        --clump-p1 5e-8 \\
        --clump-p2 1e-5 \\
        --clump-r2 0.1 \\
        --clump-kb 500 \\
        --out ${prefix} \\
        $args

    # Generate ranges for fine-mapping
    awk 'NR>1 {print \$1, \$4-500000, \$4+500000, \$3}' ${prefix}.clumped | \\
        awk '{if(\$2<0) \$2=0; print}' > ${prefix}.clumped.ranges

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | head -1 | sed 's/PLINK v//' | sed 's/ .*//')
    END_VERSIONS
    """
}
