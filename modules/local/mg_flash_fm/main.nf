process MG_FLASH_FM {
    tag "$trait_meta.trait"
    label 'process_high'

    // MG-FLASH-FM is implemented in R
    conda "conda-forge::r-base=4.3 conda-forge::r-data.table conda-forge::r-matrix"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/mg-flash-fm:latest' :
        'your-registry/mg-flash-fm:latest' }"

    input:
    tuple val(trait_meta), val(ancestries), path(sumstats_files)
    path ld_reference_dir
    val window_kb
    val max_causal

    output:
    tuple val(trait_meta), path("${prefix}.mg_flash_fm.results.tsv.gz"), emit: results
    tuple val(trait_meta), path("${prefix}.mg_flash_fm.credible_sets.tsv"), emit: credible_sets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}.multi_ancestry"
    """
    # MG-FLASH-FM: Multi-trait, multi-ancestry fine-mapping
    # Designed for related traits across multiple ancestries
    # Leverages shared genetic architecture while accounting for LD differences

    # Prepare input configuration
    cat > mg_flash_fm_config.json << EOF
    {
        "traits": ["${trait_meta.trait}"],
        "ancestries": [${ancestries.collect { '"' + it + '"' }.join(', ')}],
        "sumstats_files": [${sumstats_files.collect { '"' + it + '"' }.join(', ')}],
        "ld_reference_dir": "${ld_reference_dir}",
        "window_kb": ${window_kb},
        "max_causal": ${max_causal},
        "output_prefix": "${prefix}"
    }
    EOF

    # Run MG-FLASH-FM
    Rscript run_mg_flash_fm.R \\
        --config mg_flash_fm_config.json \\
        --threads ${task.cpus} \\
        $args

    # Compress results
    gzip ${prefix}.mg_flash_fm.results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mg_flash_fm: "1.0.0"
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
