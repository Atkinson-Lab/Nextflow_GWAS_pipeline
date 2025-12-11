process REGIONAL_PLOT {
    tag "$meta.trait - $meta.ancestry"
    label 'process_low'

    conda "conda-forge::r-base=4.3 bioconda::locuszoomr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/locuszoom:latest' :
        'your-registry/locuszoom:latest' }"

    input:
    tuple val(meta), path(sumstats)
    path gene_annotation

    output:
    tuple val(meta), path("${prefix}/*.regional.png"), emit: plots
    tuple val(meta), path("${prefix}/*.regional.pdf"), emit: pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${meta.ancestry}"
    """
    mkdir -p ${prefix}

    # Create regional association plots (LocusZoom-style)
    Rscript create_regional_plots.R \\
        --sumstats ${sumstats} \\
        --genes ${gene_annotation} \\
        --output-dir ${prefix} \\
        --ancestry ${meta.ancestry} \\
        --window-kb 500 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
