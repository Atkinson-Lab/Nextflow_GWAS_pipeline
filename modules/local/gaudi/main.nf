process GAUDI {
    tag "$trait_meta.trait"
    label 'process_high'

    // GAUDI: Local ancestry-informed PRS
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gaudi:latest' :
        'your-registry/gaudi:latest' }"

    input:
    tuple val(trait_meta), path(sumstats_files), val(sources), path(local_ancestry)
    path ld_reference_dir

    output:
    tuple val(trait_meta), path("${prefix}.gaudi.weights.txt"), emit: weights
    tuple val(trait_meta), path("${prefix}.gaudi.local_weights.txt"), emit: local_weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}"
    """
    # GAUDI: Generalized Ancestry-informed PRS with Unbiased Design
    # Uses local ancestry to weight PRS contributions
    # Particularly useful for admixed populations

    # Prepare input configuration
    cat > gaudi_config.json << EOF
    {
        "trait": "${trait_meta.trait}",
        "sumstats": [${sumstats_files.collect { '"' + it + '"' }.join(', ')}],
        "local_ancestry": "${local_ancestry}",
        "ld_reference": "${ld_reference_dir}",
        "output_prefix": "${prefix}"
    }
    EOF

    # Run GAUDI
    python run_gaudi.py \\
        --config gaudi_config.json \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gaudi: "1.0.0"
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
