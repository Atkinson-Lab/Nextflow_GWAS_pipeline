process COLOC {
    tag "$meta.trait - $qtl_type"
    label 'process_medium'

    conda "bioconda::r-coloc=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-coloc:5.2.0' :
        'quay.io/biocontainers/r-coloc:5.2.0' }"

    input:
    tuple val(meta), path(gwas_sumstats), val(qtl_type), path(qtl_sumstats)
    val p1
    val p2
    val p12

    output:
    tuple val(meta), path("${prefix}.coloc_results.tsv"), emit: results
    tuple val(meta), path("${prefix}.coloc_sensitivity.pdf"), emit: sensitivity
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${qtl_type}"
    """
    # Colocalization analysis using coloc package
    # Tests whether GWAS and QTL signals share the same causal variant

    Rscript run_coloc.R \\
        --gwas ${gwas_sumstats} \\
        --qtl ${qtl_sumstats} \\
        --qtl-type ${qtl_type} \\
        --p1 ${p1} \\
        --p2 ${p2} \\
        --p12 ${p12} \\
        --output-prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coloc: \$(Rscript -e "cat(as.character(packageVersion('coloc')))")
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
