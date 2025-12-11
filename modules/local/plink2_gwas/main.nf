process PLINK2_GWAS {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype)
    val covariate_cols
    val gwas_model

    output:
    tuple val(meta), path("${prefix}.${meta.trait}.glm.*.gz"), emit: summary_stats
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ? "--covar-name ${covariate_cols}" : ""
    def glm_type = meta.binary ? "logistic" : "linear"
    def model_arg = gwas_model == 'dominant' ? "dominant" : (gwas_model == 'recessive' ? "recessive" : "")
    """
    plink2 \\
        --bfile ${bed.baseName} \\
        --pheno ${phenotype} \\
        --pheno-name ${meta.trait} \\
        ${covar_arg} \\
        --glm ${model_arg} hide-covar \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        $args \\
        2>&1 | tee ${prefix}.log

    # Compress output
    gzip ${prefix}.${meta.trait}.glm.*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
