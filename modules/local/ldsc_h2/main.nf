process LDSC_H2 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_medium'

    conda "bioconda::ldsc=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--pyhdfd78af_2' :
        'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(munged_sumstats)
    path ld_reference_dir

    output:
    tuple val(meta), path("${prefix}.h2.log"), emit: h2
    tuple val(meta), path("${prefix}.h2_results.tsv"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def ancestry_ld = meta.ancestry ?: 'EUR'
    """
    # LDSC heritability estimation
    # Uses LD Score regression to estimate SNP-heritability

    ldsc.py \\
        --h2 ${munged_sumstats} \\
        --ref-ld-chr ${ld_reference_dir}/${ancestry_ld}/ld. \\
        --w-ld-chr ${ld_reference_dir}/${ancestry_ld}/weights. \\
        --out ${prefix}.h2 \\
        $args

    # Parse results to TSV
    python3 << 'PYEOF'
import re

with open("${prefix}.h2.log", 'r') as f:
    log = f.read()

h2 = re.search(r'Total Observed scale h2: ([0-9.]+) \\(([0-9.]+)\\)', log)
intercept = re.search(r'Intercept: ([0-9.]+) \\(([0-9.]+)\\)', log)
lambda_gc = re.search(r'Lambda GC: ([0-9.]+)', log)
mean_chi2 = re.search(r'Mean Chi\\^2: ([0-9.]+)', log)

with open("${prefix}.h2_results.tsv", 'w') as f:
    f.write("trait\\tancestry\\th2\\th2_se\\tintercept\\tintercept_se\\tlambda_gc\\tmean_chi2\\n")
    f.write(f"${meta.trait}\\t${meta.ancestry}\\t{h2.group(1) if h2 else 'NA'}\\t{h2.group(2) if h2 else 'NA'}\\t")
    f.write(f"{intercept.group(1) if intercept else 'NA'}\\t{intercept.group(2) if intercept else 'NA'}\\t")
    f.write(f"{lambda_gc.group(1) if lambda_gc else 'NA'}\\t{mean_chi2.group(1) if mean_chi2 else 'NA'}\\n")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ldsc: \$(ldsc.py --version 2>&1 | head -1 || echo "1.0.1")
    END_VERSIONS
    """
}
