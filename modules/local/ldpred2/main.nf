process LDPRED2 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    conda "bioconda::r-bigsnpr=1.12.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-bigsnpr:1.12.2' :
        'quay.io/biocontainers/r-bigsnpr:1.12.2' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_reference_dir

    output:
    tuple val(meta), path("${prefix}.ldpred2.weights.tsv"), emit: weights
    tuple val(meta), path("${prefix}.ldpred2.h2.txt"), emit: h2_estimate
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    # LDpred2: Bayesian PRS method using LD information
    # Implements auto model for automatic hyperparameter tuning

    Rscript run_ldpred2.R \\
        --sumstats ${sumstats} \\
        --ld-dir ${ld_reference_dir}/${meta.ancestry} \\
        --output-prefix ${prefix} \\
        --model auto \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ldpred2: \$(Rscript -e "cat(as.character(packageVersion('bigsnpr')))")
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
