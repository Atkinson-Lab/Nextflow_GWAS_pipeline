process MR_MEGA_ENV {
    tag "$trait_meta.trait"
    label 'process_medium'

    // MR-MEGA container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/mr-mega:latest' :
        'your-registry/mr-mega:latest' }"

    input:
    tuple val(trait_meta), val(ancestries), path(sumstats_files)
    val env_variable

    output:
    tuple val(trait_meta), path("${prefix}.mr_mega_env.txt.gz"), emit: meta_results
    tuple val(trait_meta), path("${prefix}.gxe_heterogeneity.txt"), emit: gxe_heterogeneity
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}.meta_env"
    """
    # MR-MEGA with environmental covariate
    # Extends MR-MEGA to include study-level environmental variables
    # Useful for detecting GxE interactions across ancestries

    # Create input specification file with environmental variable
    echo "study ancestry file env_${env_variable}" > mr_mega_env_input.txt
    for i in \$(seq 0 \$((${ancestries.size()} - 1))); do
        anc=\$(echo "${ancestries.join(' ')}" | cut -d' ' -f\$((i+1)))
        file=\$(echo "${sumstats_files.join(' ')}" | cut -d' ' -f\$((i+1)))
        # Environmental variable values would be extracted from study metadata
        echo "\${anc}_study \${anc} \${file} \${ENV_VALUE}" >> mr_mega_env_input.txt
    done

    # Run MR-MEGA-env
    MR-MEGA \\
        --input mr_mega_env_input.txt \\
        --output ${prefix}.mr_mega_env.txt \\
        --pc_axes 4 \\
        --env ${env_variable} \\
        --qt \\
        $args

    # Extract GxE heterogeneity statistics
    extract_mr_mega_gxe.py \\
        --input ${prefix}.mr_mega_env.txt \\
        --output ${prefix}.gxe_heterogeneity.txt

    gzip ${prefix}.mr_mega_env.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mr-mega: \$(MR-MEGA --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
