process SPA_COX {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    // SPACox for time-to-event analysis in biobank-scale data
    // Handles population structure and admixed populations
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://wzhou88/saige:1.3.0' :
        'wzhou88/saige:1.3.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype)
    val covariate_cols
    val time_col       // Time-to-event column
    val event_col      // Event indicator column (1=event, 0=censored)

    output:
    tuple val(meta), path("${prefix}.spacox.results.tsv.gz"), emit: summary_stats
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ?: ""
    """
    # SPACox: Saddlepoint Approximation for Cox proportional hazards
    # Particularly useful for:
    # - Time-to-event phenotypes (survival analysis)
    # - Admixed populations with case-control imbalance
    # - Large biobank-scale data

    # Step 1: Fit null Cox model with saddlepoint approximation
    Rscript /usr/local/bin/step1_fitNULLGLMM.R \\
        --plinkFile=${bed.baseName} \\
        --phenoFile=${phenotype} \\
        --phenoCol=${event_col} \\
        --eventTimeCol=${time_col} \\
        --covarColList=${covar_arg} \\
        --traitType=survival \\
        --outputPrefix=${prefix}.step1 \\
        --nThreads=${task.cpus} \\
        --useSparseGRMtoFitNULL=TRUE \\
        $args \\
        2>&1 | tee ${prefix}.log

    # Step 2: Run association testing
    Rscript /usr/local/bin/step2_SPAtests.R \\
        --bedFile=${bed} \\
        --bimFile=${bim} \\
        --famFile=${fam} \\
        --GMMATmodelFile=${prefix}.step1.rda \\
        --varianceRatioFile=${prefix}.step1.varianceRatio.txt \\
        --SAIGEOutputFile=${prefix}.spacox.results.tsv \\
        --minMAC=20 \\
        --LOCO=TRUE \\
        --is_Firth_beta=TRUE \\
        2>&1 | tee -a ${prefix}.log

    gzip ${prefix}.spacox.results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spacox: "via SAIGE 1.3.0"
    END_VERSIONS
    """
}

process SPA_COX_STEP1 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://wzhou88/saige:1.3.0' :
        'wzhou88/saige:1.3.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(phenotype)
    val covariate_cols
    val time_col
    val event_col

    output:
    tuple val(meta), path("${prefix}.rda"), path("${prefix}.varianceRatio.txt"), emit: model
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}.spacox_step1"
    def covar_arg = covariate_cols ?: ""
    """
    Rscript /usr/local/bin/step1_fitNULLGLMM.R \\
        --plinkFile=${bed.baseName} \\
        --phenoFile=${phenotype} \\
        --phenoCol=${event_col} \\
        --eventTimeCol=${time_col} \\
        --covarColList=${covar_arg} \\
        --sampleIDColinphenoFile=IID \\
        --traitType=survival \\
        --outputPrefix=${prefix} \\
        --nThreads=${task.cpus} \\
        --useSparseGRMtoFitNULL=TRUE \\
        --isCateVarianceRatio=TRUE \\
        $args \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saige_spacox: "1.3.0"
    END_VERSIONS
    """
}

process SPA_COX_STEP2 {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://wzhou88/saige:1.3.0' :
        'wzhou88/saige:1.3.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(model_rda), path(variance_ratio)

    output:
    tuple val(meta), path("${prefix}.spacox.tsv.gz"), emit: summary_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    Rscript /usr/local/bin/step2_SPAtests.R \\
        --bedFile=${bed} \\
        --bimFile=${bim} \\
        --famFile=${fam} \\
        --GMMATmodelFile=${model_rda} \\
        --varianceRatioFile=${variance_ratio} \\
        --SAIGEOutputFile=${prefix}.spacox.tsv \\
        --minMAC=20 \\
        --minMAF=0.0001 \\
        --LOCO=TRUE \\
        --is_Firth_beta=TRUE \\
        --is_output_moreDetails=TRUE \\
        $args

    gzip ${prefix}.spacox.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saige_spacox: "1.3.0"
    END_VERSIONS
    """
}
