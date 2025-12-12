/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FORMAT CONVERSION MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Comprehensive format conversion for GWAS tools
    Supports: PLINK1/2, VCF, GDS, BGEN formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PLINK_TO_VCF {
    tag "$meta.id - $group_label"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10 bioconda::bcftools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-plink2-bcftools:2.00a5.10--1.19' :
        'quay.io/biocontainers/mulled-v2-plink2-bcftools:2.00a5.10--1.19' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val group_label

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Convert PLINK to VCF
    plink2 \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --export vcf bgz \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        $args

    # Index VCF
    bcftools index -t ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process VCF_TO_PLINK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_idx)

    output:
    tuple val(meta), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: plink
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --vcf ${vcf} \\
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

process MERGE_PLINK_FILES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink=1.90b7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b7.2--h031d066_0' :
        'quay.io/biocontainers/plink:1.90b7.2--h031d066_0' }"

    input:
    tuple val(meta), path(beds), path(bims), path(fams)

    output:
    tuple val(meta), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: merged
    tuple val(meta), path("${prefix}.merge_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_list = beds instanceof List ? beds : [beds]
    def bim_list = bims instanceof List ? bims : [bims]
    def fam_list = fams instanceof List ? fams : [fams]
    """
    # Create merge list file
    # First file is the reference, remaining files are in the merge list
    if [ ${bed_list.size()} -gt 1 ]; then
        for i in \$(seq 1 \$((${bed_list.size()} - 1))); do
            bed_arr=(${bed_list.join(' ')})
            bim_arr=(${bim_list.join(' ')})
            fam_arr=(${fam_list.join(' ')})
            echo "\${bed_arr[\$i]%.*} \${bim_arr[\$i]} \${fam_arr[\$i]}" >> merge_list.txt
        done

        # First file as base
        first_bed="${bed_list[0]}"
        first_bim="${bim_list[0]}"
        first_fam="${fam_list[0]}"

        plink \\
            --bed \${first_bed} \\
            --bim \${first_bim} \\
            --fam \${first_fam} \\
            --merge-list merge_list.txt \\
            --make-bed \\
            --out ${prefix} \\
            --threads ${task.cpus} \\
            $args
    else
        # Only one file, just copy
        cp ${bed_list[0]} ${prefix}.bed
        cp ${bim_list[0]} ${prefix}.bim
        cp ${fam_list[0]} ${prefix}.fam
    fi

    # Create merge log
    echo "Merged ${bed_list.size()} PLINK files" > ${prefix}.merge_log.txt
    echo "Source files:" >> ${prefix}.merge_log.txt
    for f in ${bed_list.join(' ')}; do
        echo "  - \$f" >> ${prefix}.merge_log.txt
    done
    echo "Output: ${prefix}" >> ${prefix}.merge_log.txt
    echo "Samples: \$(wc -l < ${prefix}.fam)" >> ${prefix}.merge_log.txt
    echo "Variants: \$(wc -l < ${prefix}.bim)" >> ${prefix}.merge_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version | head -1 | sed 's/PLINK v//' | sed 's/ .*//')
    END_VERSIONS
    """
}

process BGEN_TO_PLINK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bgen), path(sample)

    output:
    tuple val(meta), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: plink
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --bgen ${bgen} ref-first \\
        --sample ${sample} \\
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

process PLINK_TO_BGEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${prefix}.bgen"), path("${prefix}.sample"), emit: bgen
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --export bgen-1.2 \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}

process STANDARDIZE_SUMSTATS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(meta), path(sumstats)
    val source_format  // 'regenie', 'saige', 'bolt', 'plink2', 'genesis', 'tractor', 'metal'

    output:
    tuple val(meta), path("${prefix}.standardized.tsv.gz"), emit: standardized
    tuple val(meta), path("${prefix}.sumstats_info.json"), emit: info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import json
    import sys
    import gzip

    # Column mappings for different GWAS tools
    column_maps = {
        'regenie': {
            'CHROM': 'CHR', 'GENPOS': 'POS', 'ID': 'SNP',
            'ALLELE0': 'REF', 'ALLELE1': 'ALT', 'A1FREQ': 'MAF',
            'BETA': 'BETA', 'SE': 'SE', 'LOG10P': 'LOG10P',
            'N': 'N', 'INFO': 'INFO'
        },
        'saige': {
            'CHR': 'CHR', 'POS': 'POS', 'SNPID': 'SNP', 'MarkerID': 'SNP',
            'Allele1': 'REF', 'Allele2': 'ALT', 'AF_Allele2': 'MAF',
            'BETA': 'BETA', 'SE': 'SE', 'p.value': 'P', 'Tstat': 'TSTAT',
            'N': 'N', 'imputationInfo': 'INFO'
        },
        'bolt': {
            'CHR': 'CHR', 'BP': 'POS', 'SNP': 'SNP',
            'ALLELE1': 'REF', 'ALLELE0': 'ALT', 'A1FREQ': 'MAF',
            'BETA': 'BETA', 'SE': 'SE', 'P_BOLT_LMM_INF': 'P', 'P_BOLT_LMM': 'P',
            'INFO': 'INFO'
        },
        'plink2': {
            '#CHROM': 'CHR', 'POS': 'POS', 'ID': 'SNP',
            'REF': 'REF', 'ALT': 'ALT', 'A1_FREQ': 'MAF',
            'BETA': 'BETA', 'SE': 'SE', 'P': 'P', 'LOG10_P': 'LOG10P',
            'OBS_CT': 'N', 'T_STAT': 'TSTAT'
        },
        'genesis': {
            'chr': 'CHR', 'pos': 'POS', 'variant.id': 'SNP',
            'ref': 'REF', 'alt': 'ALT', 'freq': 'MAF',
            'Est': 'BETA', 'Est.SE': 'SE', 'Score.pval': 'P',
            'n.obs': 'N'
        },
        'tractor': {
            'CHR': 'CHR', 'POS': 'POS', 'SNP': 'SNP',
            'REF': 'REF', 'ALT': 'ALT',
            'P_joint': 'P', 'P_het': 'P_HET'
        },
        'metal': {
            'Chromosome': 'CHR', 'Position': 'POS', 'MarkerName': 'SNP',
            'Allele1': 'REF', 'Allele2': 'ALT', 'Freq1': 'MAF',
            'Effect': 'BETA', 'StdErr': 'SE', 'P-value': 'P', 'Pvalue': 'P',
            'HetISq': 'I2', 'HetPVal': 'P_HET', 'Direction': 'DIRECTION'
        }
    }

    # Read input file
    if "${sumstats}".endswith('.gz'):
        df = pd.read_csv("${sumstats}", sep='\\t', compression='gzip', low_memory=False)
    else:
        df = pd.read_csv("${sumstats}", sep='\\t', low_memory=False)

    # Get column mapping for source format
    source = "${source_format}".lower()
    col_map = column_maps.get(source, {})

    # Standardize column names
    rename_dict = {}
    for old_col, new_col in col_map.items():
        if old_col in df.columns and new_col not in rename_dict.values():
            rename_dict[old_col] = new_col

    df = df.rename(columns=rename_dict)

    # Convert LOG10P to P if needed
    if 'LOG10P' in df.columns and 'P' not in df.columns:
        df['P'] = 10 ** (-df['LOG10P'])

    # Ensure required columns exist
    required = ['CHR', 'POS', 'SNP', 'BETA', 'SE', 'P']
    missing = [c for c in required if c not in df.columns]

    # Create info dictionary
    info = {
        'source_format': source,
        'n_variants': len(df),
        'columns_original': list(df.columns),
        'columns_standardized': list(df.columns),
        'missing_required': missing,
        'ancestry': '${meta.ancestry}' if '${meta.ancestry}' else None,
        'trait': '${meta.trait}' if hasattr(${meta}, 'trait') else None
    }

    if 'P' in df.columns:
        info['min_p'] = float(df['P'].min())
        info['n_significant_5e8'] = int((df['P'] < 5e-8).sum())
        info['n_suggestive_1e5'] = int((df['P'] < 1e-5).sum())

    # Save standardized summary stats
    df.to_csv("${prefix}.standardized.tsv.gz", sep='\\t', index=False, compression='gzip')

    # Save info
    with open("${prefix}.sumstats_info.json", 'w') as f:
        json.dump(info, f, indent=2)

    print(f"Standardized {len(df)} variants from {source} format")
    if missing:
        print(f"Warning: Missing required columns: {missing}", file=sys.stderr)

    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pd.__version__}"\\n')
    """
}

process MERGE_LOCAL_ANCESTRY {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(meta), path(msp_files)

    output:
    tuple val(meta), path("${prefix}.merged.msp.tsv"), emit: merged_msp
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def msp_list = msp_files instanceof List ? msp_files : [msp_files]
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    msp_files = "${msp_list.join(' ')}".split()

    # Read and merge MSP files
    dfs = []
    for f in msp_files:
        # MSP format: #chrom, spos, epos, sgpos, egpos, n_snps, sample1.0, sample1.1, sample2.0, ...
        df = pd.read_csv(f, sep='\\t', comment='#')
        dfs.append(df)

    if len(dfs) > 1:
        # Concatenate chromosome-wise
        merged = pd.concat(dfs, ignore_index=True)
        merged = merged.sort_values(['#chrom', 'spos'])
    else:
        merged = dfs[0]

    # Save merged file
    merged.to_csv("${prefix}.merged.msp.tsv", sep='\\t', index=False)

    print(f"Merged {len(msp_files)} MSP files: {len(merged)} segments")

    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
    """
}
