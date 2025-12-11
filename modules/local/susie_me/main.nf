process SUSIE_ME {
    tag "$trait_meta.trait"
    label 'process_high'

    // SuSIE-ME (Multi-Ethnic SuSIE)
    conda "conda-forge::r-base=4.3 bioconda::r-susie"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/susie-me:latest' :
        'your-registry/susie-me:latest' }"

    input:
    tuple val(trait_meta), val(ancestries), path(sumstats_files)
    path ld_reference_dir
    val max_causal

    output:
    tuple val(trait_meta), path("${prefix}.susie_me.results.tsv.gz"), emit: results
    tuple val(trait_meta), path("${prefix}.susie_me.credible_sets.tsv"), emit: credible_sets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}.multi_ancestry"
    """
    # SuSIE-ME: Multi-ethnic extension of SuSIE
    # For single traits or unrelated traits across ancestries
    # Uses ancestry-specific LD to improve fine-mapping

    # Prepare input files
    mkdir -p input_data

    # Create sumstats list
    echo "${sumstats_files.join('\n')}" > input_data/sumstats_list.txt
    echo "${ancestries.join('\n')}" > input_data/ancestry_list.txt

    # Run SuSIE-ME
    Rscript run_susie_me.R \\
        --sumstats-list input_data/sumstats_list.txt \\
        --ancestry-list input_data/ancestry_list.txt \\
        --ld-dir ${ld_reference_dir} \\
        --max-causal ${max_causal} \\
        --output-prefix ${prefix} \\
        --threads ${task.cpus} \\
        $args

    # Compress results
    gzip ${prefix}.susie_me.results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        susie: \$(Rscript -e "cat(as.character(packageVersion('susieR')))")
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
