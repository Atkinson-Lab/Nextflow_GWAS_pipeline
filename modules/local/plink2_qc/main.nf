process PLINK2_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val maf_threshold
    val hwe_threshold
    val geno_missing
    val mind_missing

    output:
    tuple val(meta), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: genotypes
    tuple val(meta), path("${prefix}.qc_stats.txt"), emit: qc_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.qc"
    """
    # Initial sample count
    wc -l ${fam} | cut -d' ' -f1 > initial_samples.txt
    wc -l ${bim} | cut -d' ' -f1 > initial_variants.txt

    # Run QC
    plink2 \\
        --bfile ${bed.baseName} \\
        --maf ${maf_threshold} \\
        --hwe ${hwe_threshold} \\
        --geno ${geno_missing} \\
        --mind ${mind_missing} \\
        --make-bed \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        $args

    # Final counts
    wc -l ${prefix}.fam | cut -d' ' -f1 > final_samples.txt
    wc -l ${prefix}.bim | cut -d' ' -f1 > final_variants.txt

    # Generate QC stats
    echo "QC Statistics for ${meta.id}" > ${prefix}.qc_stats.txt
    echo "================================" >> ${prefix}.qc_stats.txt
    echo "Initial samples: \$(cat initial_samples.txt)" >> ${prefix}.qc_stats.txt
    echo "Final samples: \$(cat final_samples.txt)" >> ${prefix}.qc_stats.txt
    echo "Initial variants: \$(cat initial_variants.txt)" >> ${prefix}.qc_stats.txt
    echo "Final variants: \$(cat final_variants.txt)" >> ${prefix}.qc_stats.txt
    echo "MAF threshold: ${maf_threshold}" >> ${prefix}.qc_stats.txt
    echo "HWE threshold: ${hwe_threshold}" >> ${prefix}.qc_stats.txt
    echo "Genotype missingness: ${geno_missing}" >> ${prefix}.qc_stats.txt
    echo "Sample missingness: ${mind_missing}" >> ${prefix}.qc_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
