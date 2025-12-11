process MR_MEGA {
    tag "$trait_meta.trait"
    label 'process_medium'

    // MR-MEGA container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/mr-mega:latest' :
        'your-registry/mr-mega:latest' }"

    input:
    tuple val(trait_meta), val(ancestries), path(sumstats_files)

    output:
    tuple val(trait_meta), path("${prefix}.mr_mega.txt.gz"), emit: meta_results
    tuple val(trait_meta), path("${prefix}.heterogeneity.txt"), emit: heterogeneity
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}.meta"
    """
    # Prepare input file list for MR-MEGA
    # MR-MEGA uses ancestry-specific summary statistics to model
    # allelic effects as a function of genetic distance

    # Create input specification file
    echo "study ancestry file" > mr_mega_input.txt
    for i in \$(seq 0 \$((${ancestries.size()} - 1))); do
        anc=\$(echo "${ancestries.join(' ')}" | cut -d' ' -f\$((i+1)))
        file=\$(echo "${sumstats_files.join(' ')}" | cut -d' ' -f\$((i+1)))
        echo "\${anc}_study \${anc} \${file}" >> mr_mega_input.txt
    done

    # Run MR-MEGA
    # MR-MEGA performs meta-regression of genetic association data
    # accounting for ancestry differences
    MR-MEGA \\
        --input mr_mega_input.txt \\
        --output ${prefix}.mr_mega.txt \\
        --pc_axes 4 \\
        --qt \\
        $args

    # Extract heterogeneity statistics
    # Includes ancestry-correlated and ancestry-uncorrelated heterogeneity
    extract_mr_mega_het.py \\
        --input ${prefix}.mr_mega.txt \\
        --output ${prefix}.heterogeneity.txt

    gzip ${prefix}.mr_mega.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mr-mega: \$(MR-MEGA --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
