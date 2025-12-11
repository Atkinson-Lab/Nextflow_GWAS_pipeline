process GRAF_ANCESTRY {
    tag "$meta.id"
    label 'process_medium'

    // GRAF-ANC container - you may need to build this
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/graf-anc:latest' :
        'your-registry/graf-anc:latest' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path graf_reference

    output:
    tuple val(meta), path("${prefix}.ancestry.txt"), emit: ancestry_calls
    tuple val(meta), path("${prefix}.ancestry_pcs.txt"), emit: pcs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # GRAF-ANC ancestry inference
    # Outputs 9 ancestry groups:
    # EUR (European), AFR (African), EAS (East Asian), SAS (South Asian),
    # AMR (American), MID (Middle Eastern), AAC (African American),
    # AHI (American Hispanic), HET (Heterogeneous/other)

    # Run GRAF-pop for ancestry inference
    graf -geno ${bed} \\
        -ref ${graf_reference} \\
        -out ${prefix} \\
        -threads ${task.cpus} \\
        $args

    # Process GRAF output to standardized format
    # Format: SAMPLE_ID, ANCESTRY, EUR_PROB, AFR_PROB, EAS_PROB, SAS_PROB, AMR_PROB, MID_PROB
    process_graf_output.py \\
        --input ${prefix}.txt \\
        --output ${prefix}.ancestry.txt \\
        --pcs ${prefix}.ancestry_pcs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graf: \$(graf --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}
