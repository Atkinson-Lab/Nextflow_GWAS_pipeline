process PRS_COMPARISON {
    tag "$trait_meta.trait"
    label 'process_medium'

    conda "conda-forge::r-base=4.3 conda-forge::r-data.table conda-forge::r-ggplot2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-validation:latest' :
        'your-registry/prs-validation:latest' }"

    input:
    tuple val(trait_meta), val(methods), val(ancestries), path(prs_scores)
    path validation_results

    output:
    tuple val(trait_meta), path("${prefix}.method_comparison.tsv"), emit: comparison
    tuple val(trait_meta), path("${prefix}.comparison_plot.pdf"), emit: plot
    tuple val(trait_meta), path("${prefix}.concordance.tsv"), emit: concordance
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}"
    """
    # Compare PRS methods across ancestries
    # Determines best method overall and per ancestry group
    # Calculates concordance between methods

    Rscript compare_prs_methods.R \\
        --prs-files "${prs_scores.join(',')}" \\
        --methods "${methods.join(',')}" \\
        --ancestries "${ancestries.join(',')}" \\
        --validation-dir "${validation_results ?: '.'}" \\
        --output-prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
