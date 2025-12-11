process PRS_CSX {
    tag "$trait_meta.trait"
    label 'process_high'

    conda "conda-forge::python=3.11 conda-forge::scipy conda-forge::h5py"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-csx:latest' :
        'your-registry/prs-csx:latest' }"

    input:
    tuple val(trait_meta), path(sumstats_files), val(ancestries)
    path prscsx_reference

    output:
    tuple val(trait_meta), path("${prefix}/*.prs_csx.weights.txt"), emit: weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${trait_meta.trait}"
    """
    mkdir -p ${prefix}

    # PRS-CSx: Cross-population PRS using coupled continuous shrinkage priors
    # Jointly models GWAS summary statistics from multiple populations
    # Leverages shared genetic effects while accounting for population differences

    # Prepare population-specific inputs
    pop_args=""
    for i in \$(seq 0 \$((${ancestries.size()} - 1))); do
        anc=\$(echo "${ancestries.join(' ')}" | cut -d' ' -f\$((i+1)))
        file=\$(echo "${sumstats_files.join(' ')}" | cut -d' ' -f\$((i+1)))
        pop_args="\${pop_args} --sst_file_\${anc}=\${file}"

        # Get sample size from sumstats header or metadata
        n=\$(zcat \$file | head -1000 | grep -oP 'N=\\K[0-9]+' | head -1 || echo "10000")
        pop_args="\${pop_args} --n_gwas_\${anc}=\${n}"
    done

    # Run PRS-CSx
    python PRScsx.py \\
        --ref_dir=${prscsx_reference} \\
        --bim_prefix=\${BIM_PREFIX:-NA} \\
        \${pop_args} \\
        --pop=${ancestries.join(',')} \\
        --out_dir=${prefix} \\
        --out_name=${trait_meta.trait} \\
        --phi=1e-2 \\
        --n_iter=1000 \\
        --n_burnin=500 \\
        --thin=5 \\
        --seed=42 \\
        $args

    # Rename outputs for clarity
    for f in ${prefix}/*_pst_eff_*.txt; do
        mv \$f \${f/_pst_eff_/.prs_csx.weights.}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prs-csx: "1.0.0"
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
