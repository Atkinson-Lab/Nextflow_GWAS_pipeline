process SAIGE_STEP2 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_medium'

    conda "bioconda::r-saige=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-saige:1.3.0' :
        'wzhou88/saige:1.3.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype), path(model_rda), path(variance_ratio)
    val gwas_model

    output:
    tuple val(meta), path("${prefix}.saige.txt.gz"), emit: summary_stats
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    Rscript /usr/local/bin/step2_SPAtests.R \\
        --bedFile=${bed} \\
        --bimFile=${bim} \\
        --famFile=${fam} \\
        --GMMATmodelFile=${model_rda} \\
        --varianceRatioFile=${variance_ratio} \\
        --SAIGEOutputFile=${prefix}.saige.txt \\
        --minMAC=20 \\
        --minMAF=0.0001 \\
        --LOCO=TRUE \\
        $args \\
        2>&1 | tee ${prefix}.log

    gzip ${prefix}.saige.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saige: \$(Rscript -e "library(SAIGE); cat(as.character(packageVersion('SAIGE')))")
    END_VERSIONS
    """
}
