process REGENIE_STEP2 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_medium'

    conda "bioconda::regenie=3.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/regenie:3.4.1--h2810c46_0' :
        'quay.io/biocontainers/regenie:3.4.1--h2810c46_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype), path(predictions)
    val covariate_cols
    val gwas_model

    output:
    tuple val(meta), path("${prefix}_${meta.trait}.regenie.gz"), emit: summary_stats
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--firth --approx --firth-se --minMAC 20'
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ? "--covarColList ${covariate_cols}" : ""
    def binary_arg = meta.binary ? "--bt" : ""
    def model_arg = gwas_model == 'dominant' ? "--dominance" : (gwas_model == 'recessive' ? "--recessive" : "")
    def pred_list = predictions.find { it.name.endsWith('_pred.list') }
    """
    regenie \\
        --step 2 \\
        --bed ${bed.baseName} \\
        --phenoFile ${phenotype} \\
        --phenoCol ${meta.trait} \\
        ${covar_arg} \\
        ${binary_arg} \\
        --pred ${pred_list} \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${model_arg} \\
        $args

    # Compress output
    gzip -c ${prefix}_${meta.trait}.regenie > ${prefix}_${meta.trait}.regenie.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version | sed 's/.*v//')
    END_VERSIONS
    """
}
