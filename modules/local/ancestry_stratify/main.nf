process ANCESTRY_STRATIFY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(ancestry_calls)

    output:
    tuple val(meta_out), path("*.bed"), path("*.bim"), path("*.fam"), emit: stratified_genotypes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // Output will be multiple files, one per ancestry group
    // meta_out will have ancestry field added
    """
    # Get unique ancestry groups
    cut -f2 ${ancestry_calls} | tail -n+2 | sort -u > ancestry_groups.txt

    # Stratify by each ancestry group
    while read ancestry; do
        # Extract sample IDs for this ancestry
        awk -v anc="\$ancestry" '\$2 == anc {print \$1, \$1}' ${ancestry_calls} > \${ancestry}_samples.txt

        # Skip if no samples
        if [ ! -s \${ancestry}_samples.txt ]; then
            continue
        fi

        # Extract samples using PLINK2
        plink2 \\
            --bfile ${bed.baseName} \\
            --keep \${ancestry}_samples.txt \\
            --make-bed \\
            --out ${prefix}.\${ancestry} \\
            --threads ${task.cpus} \\
            $args

    done < ancestry_groups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """

    // Note: The meta_out with ancestry field is handled in the subworkflow
    // by parsing the output filenames
}
