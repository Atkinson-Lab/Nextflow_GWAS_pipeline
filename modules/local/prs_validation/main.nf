process PRS_VALIDATION {
    tag "$meta.id - $meta.method"
    label 'process_medium'

    conda "conda-forge::r-base=4.3 conda-forge::r-data.table conda-forge::r-ggplot2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-validation:latest' :
        'your-registry/prs-validation:latest' }"

    input:
    tuple val(meta), path(prs_scores)
    path validation_cohort

    output:
    tuple val(meta), path("${prefix}.validation_results.tsv"), emit: results
    tuple val(meta), path("${prefix}.validation_metrics.json"), emit: metrics
    tuple val(meta), path("${prefix}.validation_plot.pdf"), emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${meta.method}.${meta.ancestry ?: 'multi'}"
    """
    # Validate PRS performance in target cohort
    # Calculates R², AUC (for binary), concordance metrics

    Rscript validate_prs.R \\
        --prs-scores ${prs_scores} \\
        --phenotype ${validation_cohort} \\
        --trait ${meta.trait} \\
        --output-prefix ${prefix} \\
        --binary ${meta.binary ?: 'false'} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
