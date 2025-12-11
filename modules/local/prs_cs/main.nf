process PRS_CS {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    conda "conda-forge::python=3.11 conda-forge::scipy conda-forge::h5py"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-cs:latest' :
        'your-registry/prs-cs:latest' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_reference_dir

    output:
    tuple val(meta), path("${prefix}.prs_cs.weights.txt"), emit: weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    # PRS-CS: Polygenic prediction with continuous shrinkage priors
    # Single-ancestry PRS method with Bayesian regression

    # Get sample size
    n=\$(zcat ${sumstats} | head -1000 | grep -oP 'N=\\K[0-9]+' | head -1 || echo "10000")

    python PRScs.py \\
        --ref_dir=${ld_reference_dir}/${meta.ancestry} \\
        --bim_prefix=\${BIM_PREFIX:-NA} \\
        --sst_file=${sumstats} \\
        --n_gwas=\${n} \\
        --out_dir=. \\
        --out_name=${prefix} \\
        --phi=1e-2 \\
        --n_iter=1000 \\
        --n_burnin=500 \\
        --thin=5 \\
        --seed=42 \\
        $args

    # Rename output
    mv ${prefix}_pst_eff_*.txt ${prefix}.prs_cs.weights.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prs-cs: "1.0.0"
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
