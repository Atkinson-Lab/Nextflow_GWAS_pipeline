process SOFTWARE_VERSIONS {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path versions

    output:
    path "software_versions.yml", emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${versions} > software_versions.yml

    # Create MultiQC-compatible versions file
    echo "id: 'software_versions'" > software_versions_mqc.yml
    echo "section_name: 'Software Versions'" >> software_versions_mqc.yml
    echo "plot_type: 'html'" >> software_versions_mqc.yml
    echo "data: |" >> software_versions_mqc.yml
    echo "  <dl class='dl-horizontal'>" >> software_versions_mqc.yml
    cat ${versions} | while read line; do
        echo "    <dt>\$line</dt>" >> software_versions_mqc.yml
    done
    echo "  </dl>" >> software_versions_mqc.yml
    """
}
