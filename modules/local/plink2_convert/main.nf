process PLINK2_CONVERT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(genotype_files)

    output:
    tuple val(meta), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: plink1
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def input_cmd = ""
    if (meta.format == 'plink2') {
        input_cmd = "--pfile ${genotype_files[0].baseName}"
    } else if (meta.format == 'bgen') {
        input_cmd = "--bgen ${genotype_files[0]} --sample ${genotype_files[1]}"
    } else if (meta.format == 'vcf') {
        input_cmd = "--vcf ${genotype_files[0]}"
    }

    """
    plink2 \\
        $input_cmd \\
        --make-bed \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
