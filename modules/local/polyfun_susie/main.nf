process POLYFUN_SUSIE {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    conda "conda-forge::python=3.11 conda-forge::numpy conda-forge::scipy conda-forge::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/polyfun:latest' :
        'your-registry/polyfun:latest' }"

    input:
    tuple val(meta), path(munged_sumstats)
    path ld_reference_dir
    val max_causal

    output:
    tuple val(meta), path("${prefix}/*.polyfun_susie.tsv.gz"), emit: results
    tuple val(meta), path("${prefix}/*.credible_sets.tsv"), emit: credible_sets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    mkdir -p ${prefix}

    # PolyFun+SuSIE fine-mapping
    # Uses functionally-informed prior probabilities from PolyFun
    # Combined with SuSIE for fine-mapping with multiple causal variants

    # Run PolyFun to compute prior probabilities
    python -m polyfun.polyfun \\
        --sumstats ${munged_sumstats} \\
        --ref-ld-chr ${ld_reference_dir}/${meta.ancestry}/ld. \\
        --w-ld-chr ${ld_reference_dir}/${meta.ancestry}/weights. \\
        --output-prefix ${prefix}/polyfun \\
        $args

    # Run SuSIE with PolyFun priors
    python -m polyfun.finemapper \\
        --sumstats ${prefix}/polyfun.sumstats.parquet \\
        --method susie \\
        --max-num-causal ${max_causal} \\
        --n ${meta.n_samples ?: 10000} \\
        --ld ${ld_reference_dir}/${meta.ancestry} \\
        --out-prefix ${prefix}/${meta.trait} \\
        --allow-missing \\
        --finemap-all

    # Extract credible sets
    for f in ${prefix}/*.susie.tsv; do
        extract_credible_sets.py \\
            --input \$f \\
            --output \${f%.susie.tsv}.credible_sets.tsv \\
            --pip-threshold 0.95
    done

    # Compress results
    gzip ${prefix}/*.susie.tsv
    mv ${prefix}/*.susie.tsv.gz ${prefix}/*.polyfun_susie.tsv.gz || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polyfun: \$(python -c "import polyfun; print(polyfun.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
