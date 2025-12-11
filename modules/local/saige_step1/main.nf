process SAIGE_STEP1 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    conda "bioconda::r-saige=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-saige:1.3.0' :
        'wzhou88/saige:1.3.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype)
    val covariate_cols

    output:
    tuple val(meta), path("${prefix}.rda"), path("${prefix}.varianceRatio.txt"), emit: model
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}.step1"
    def covar_arg = covariate_cols ? "covarColList=${covariate_cols}" : ""
    def trait_type = meta.binary ? "binary" : "quantitative"
    """
    Rscript /usr/local/bin/step1_fitNULLGLMM.R \\
        --plinkFile=${bed.baseName} \\
        --phenoFile=${phenotype} \\
        --phenoCol=${meta.trait} \\
        ${covar_arg} \\
        --traitType=${trait_type} \\
        --outputPrefix=${prefix} \\
        --nThreads=${task.cpus} \\
        $args \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saige: \$(Rscript -e "library(SAIGE); cat(as.character(packageVersion('SAIGE')))")
    END_VERSIONS
    """
}
