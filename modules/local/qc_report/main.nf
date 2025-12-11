process QC_REPORT {
    tag "qc_report"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    path qc_stats_files

    output:
    path "qc_summary_report.html", emit: report
    path "qc_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_qc_report.py \\
        --input ${qc_stats_files.join(' ')} \\
        --output qc_summary_report.html \\
        --summary qc_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
