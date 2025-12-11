process REGENIE_STEP1 {
    tag "$meta.id - $meta.ancestry"
    label 'process_high'

    conda "bioconda::regenie=3.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/regenie:3.4.1--h2810c46_0' :
        'quay.io/biocontainers/regenie:3.4.1--h2810c46_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype)
    val covariate_cols

    output:
    tuple val(meta), path("${prefix}_pred.list"), path("${prefix}_*.loco.gz"), emit: predictions
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--loocv --bsize 1000 --lowmem'
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.step1"
    def covar_arg = covariate_cols ? "--covarColList ${covariate_cols}" : ""
    """
    regenie \\
        --step 1 \\
        --bed ${bed.baseName} \\
        --phenoFile ${phenotype} \\
        ${covar_arg} \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version | sed 's/.*v//')
    END_VERSIONS
    """
}
