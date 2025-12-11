process LDSC_RG {
    tag "$meta.trait - $meta.anc1 vs $meta.anc2"
    label 'process_medium'

    conda "bioconda::ldsc=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--pyhdfd78af_2' :
        'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(munged1), path(munged2)
    path ld_reference_dir

    output:
    tuple val(meta), path("${prefix}.rg.log"), emit: rg
    tuple val(meta), path("${prefix}.rg_results.tsv"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${meta.anc1}_vs_${meta.anc2}"
    """
    # Cross-ancestry genetic correlation using LDSC
    # Note: For cross-ancestry rg, need appropriate cross-ancestry LD scores

    ldsc.py \\
        --rg ${munged1},${munged2} \\
        --ref-ld-chr ${ld_reference_dir}/EUR/ld. \\
        --w-ld-chr ${ld_reference_dir}/EUR/weights. \\
        --out ${prefix}.rg \\
        $args

    # Parse results
    python3 << 'PYEOF'
import re

with open("${prefix}.rg.log", 'r') as f:
    log = f.read()

rg = re.search(r'Genetic Correlation: ([0-9.-]+) \\(([0-9.]+)\\)', log)
p = re.search(r'P: ([0-9.e-]+)', log)

with open("${prefix}.rg_results.tsv", 'w') as f:
    f.write("trait\\tancestry1\\tancestry2\\trg\\trg_se\\tp\\n")
    f.write(f"${meta.trait}\\t${meta.anc1}\\t${meta.anc2}\\t")
    f.write(f"{rg.group(1) if rg else 'NA'}\\t{rg.group(2) if rg else 'NA'}\\t")
    f.write(f"{p.group(1) if p else 'NA'}\\n")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ldsc: \$(ldsc.py --version 2>&1 | head -1 || echo "1.0.1")
    END_VERSIONS
    """
}
