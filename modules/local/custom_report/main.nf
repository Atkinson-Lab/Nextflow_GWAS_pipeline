process CUSTOM_REPORT {
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::jinja2 conda-forge::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/report-generator:latest' :
        'your-registry/report-generator:latest' }"

    input:
    path gwas_results
    path meta_results
    path fm_results
    path h2_results
    path prs_results

    output:
    path "ancestry_gwas_report.html", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_pipeline_report.py \\
        --gwas-results ${gwas_results} \\
        --meta-results ${meta_results ?: 'NA'} \\
        --fm-results ${fm_results ?: 'NA'} \\
        --h2-results ${h2_results ?: 'NA'} \\
        --prs-results ${prs_results ?: 'NA'} \\
        --output ancestry_gwas_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
