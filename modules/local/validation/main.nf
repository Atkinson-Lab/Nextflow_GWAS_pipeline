/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input validation and decision-making checks for the GWAS pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VALIDATE_GWAS_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path phenotype

    output:
    tuple val(meta), path("${prefix}.validation_report.json"), emit: report
    tuple val(meta), path("${prefix}.sample_overlap.txt"), emit: sample_overlap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry ?: 'all'}"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import json
    import sys
    import os

    validation = {
        'status': 'PASS',
        'warnings': [],
        'errors': [],
        'meta': {
            'id': '${meta.id}',
            'ancestry': '${meta.ancestry}' if '${meta.ancestry}' else None
        }
    }

    # =========================================================================
    # VALIDATE PLINK FILES
    # =========================================================================

    # Check BED file exists and has size
    bed_size = os.path.getsize("${bed}")
    if bed_size == 0:
        validation['errors'].append("BED file is empty")
        validation['status'] = 'FAIL'
    validation['bed_size_mb'] = round(bed_size / 1024 / 1024, 2)

    # Read BIM file (variant info)
    try:
        bim = pd.read_csv("${bim}", sep='\\t', header=None,
                         names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
        validation['n_variants'] = len(bim)
        validation['chromosomes'] = sorted(bim['chr'].unique().tolist())

        # Check for duplicate variants
        dup_snps = bim['snp'].duplicated().sum()
        if dup_snps > 0:
            validation['warnings'].append(f"{dup_snps} duplicate SNP IDs found")
        validation['duplicate_snps'] = int(dup_snps)

        # Check for missing positions
        missing_pos = (bim['pos'] == 0).sum()
        if missing_pos > 0:
            validation['warnings'].append(f"{missing_pos} variants with position 0")

    except Exception as e:
        validation['errors'].append(f"Failed to read BIM file: {str(e)}")
        validation['status'] = 'FAIL'

    # Read FAM file (sample info)
    try:
        fam = pd.read_csv("${fam}", sep='\\s+', header=None,
                         names=['fid', 'iid', 'pid', 'mid', 'sex', 'pheno'])
        validation['n_samples'] = len(fam)
        geno_samples = set(fam['iid'].astype(str))

        # Check for duplicate samples
        dup_samples = fam['iid'].duplicated().sum()
        if dup_samples > 0:
            validation['errors'].append(f"{dup_samples} duplicate sample IDs")
            validation['status'] = 'FAIL'

        # Check sex coding
        sex_values = fam['sex'].unique()
        validation['sex_values'] = sex_values.tolist()
        if not set(sex_values).issubset({0, 1, 2}):
            validation['warnings'].append(f"Non-standard sex coding: {sex_values}")

    except Exception as e:
        validation['errors'].append(f"Failed to read FAM file: {str(e)}")
        validation['status'] = 'FAIL'
        geno_samples = set()

    # =========================================================================
    # VALIDATE PHENOTYPE FILE
    # =========================================================================

    try:
        # Auto-detect separator
        with open("${phenotype}", 'r') as f:
            first_line = f.readline()
        sep = '\\t' if '\\t' in first_line else ','

        pheno = pd.read_csv("${phenotype}", sep=sep)
        validation['n_pheno_samples'] = len(pheno)
        validation['pheno_columns'] = pheno.columns.tolist()

        # Find sample ID column
        id_cols = ['IID', 'iid', 'sample_id', 'SAMPLE_ID', 'ID', 'id', 'SampleID']
        sample_col = None
        for col in id_cols:
            if col in pheno.columns:
                sample_col = col
                break

        if sample_col is None:
            # Try first column
            sample_col = pheno.columns[0]
            validation['warnings'].append(f"Using first column as sample ID: {sample_col}")

        pheno_samples = set(pheno[sample_col].astype(str))
        validation['pheno_sample_col'] = sample_col

        # Check sample overlap
        overlap = geno_samples & pheno_samples
        geno_only = geno_samples - pheno_samples
        pheno_only = pheno_samples - geno_samples

        validation['sample_overlap'] = {
            'n_overlap': len(overlap),
            'n_geno_only': len(geno_only),
            'n_pheno_only': len(pheno_only),
            'overlap_pct': round(100 * len(overlap) / len(geno_samples), 2) if geno_samples else 0
        }

        if len(overlap) == 0:
            validation['errors'].append("No sample overlap between genotype and phenotype files!")
            validation['status'] = 'FAIL'
        elif len(overlap) < len(geno_samples) * 0.5:
            validation['warnings'].append(f"Less than 50% sample overlap ({len(overlap)}/{len(geno_samples)})")

        # Save overlapping samples
        with open("${prefix}.sample_overlap.txt", 'w') as f:
            for s in sorted(overlap):
                f.write(f"{s}\\n")

        # Check for phenotype columns with high missingness
        numeric_cols = pheno.select_dtypes(include=['number']).columns
        for col in numeric_cols:
            miss_pct = pheno[col].isna().mean() * 100
            if miss_pct > 50:
                validation['warnings'].append(f"Column '{col}' has {miss_pct:.1f}% missing values")

    except Exception as e:
        validation['errors'].append(f"Failed to read phenotype file: {str(e)}")
        validation['status'] = 'FAIL'
        with open("${prefix}.sample_overlap.txt", 'w') as f:
            pass  # Create empty file

    # =========================================================================
    # FINAL STATUS
    # =========================================================================

    if validation['errors']:
        validation['status'] = 'FAIL'
    elif validation['warnings']:
        validation['status'] = 'WARN'

    # Save validation report
    with open("${prefix}.validation_report.json", 'w') as f:
        json.dump(validation, f, indent=2, default=str)

    # Print summary
    print(f"Validation Status: {validation['status']}")
    print(f"  Variants: {validation.get('n_variants', 'N/A')}")
    print(f"  Samples (genotype): {validation.get('n_samples', 'N/A')}")
    print(f"  Samples (phenotype): {validation.get('n_pheno_samples', 'N/A')}")
    if 'sample_overlap' in validation:
        print(f"  Sample overlap: {validation['sample_overlap']['n_overlap']} ({validation['sample_overlap']['overlap_pct']}%)")

    if validation['warnings']:
        print("  Warnings:")
        for w in validation['warnings']:
            print(f"    - {w}")

    if validation['errors']:
        print("  Errors:", file=sys.stderr)
        for e in validation['errors']:
            print(f"    - {e}", file=sys.stderr)
        sys.exit(1)

    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
    """
}

process VALIDATE_ANCESTRY_CALLS {
    tag "ancestry_validation"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    path ancestry_file
    path fam_file
    val expected_groups

    output:
    path "ancestry_validation.json", emit: report
    path "ancestry_counts.tsv", emit: counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import json
    import sys

    expected = "${expected_groups}".split(',')

    validation = {
        'status': 'PASS',
        'warnings': [],
        'errors': [],
        'expected_groups': expected
    }

    # Read ancestry calls
    try:
        anc = pd.read_csv("${ancestry_file}", sep='\\t')

        # Find columns
        id_cols = ['IID', 'iid', 'sample_id', 'SAMPLE_ID', 'ID']
        anc_cols = ['ancestry', 'ANCESTRY', 'predicted_ancestry', 'population', 'POP']

        id_col = None
        for col in id_cols:
            if col in anc.columns:
                id_col = col
                break

        ancestry_col = None
        for col in anc_cols:
            if col in anc.columns:
                ancestry_col = col
                break

        if id_col is None or ancestry_col is None:
            validation['errors'].append(f"Could not identify ID or ancestry columns. Found: {anc.columns.tolist()}")
            validation['status'] = 'FAIL'
        else:
            # Count ancestry groups
            counts = anc[ancestry_col].value_counts()
            validation['ancestry_counts'] = counts.to_dict()
            validation['n_samples'] = len(anc)
            validation['unique_ancestries'] = counts.index.tolist()

            # Check for expected groups
            found = set(counts.index)
            expected_set = set(expected)
            missing_expected = expected_set - found
            unexpected = found - expected_set

            if missing_expected:
                validation['warnings'].append(f"Expected ancestries not found: {missing_expected}")
            if unexpected:
                validation['warnings'].append(f"Unexpected ancestries found: {unexpected}")

            # Check minimum sample size
            for anc_group, count in counts.items():
                if count < 100:
                    validation['warnings'].append(f"Ancestry '{anc_group}' has only {count} samples (min recommended: 100)")

            # Save counts
            counts_df = pd.DataFrame({'ancestry': counts.index, 'count': counts.values})
            counts_df.to_csv("ancestry_counts.tsv", sep='\\t', index=False)

    except Exception as e:
        validation['errors'].append(f"Failed to read ancestry file: {str(e)}")
        validation['status'] = 'FAIL'
        pd.DataFrame().to_csv("ancestry_counts.tsv", sep='\\t', index=False)

    # Read FAM to check sample overlap
    try:
        fam = pd.read_csv("${fam_file}", sep='\\s+', header=None,
                         names=['fid', 'iid', 'pid', 'mid', 'sex', 'pheno'])
        fam_samples = set(fam['iid'].astype(str))

        if 'id_col' in dir() and id_col:
            anc_samples = set(anc[id_col].astype(str))
            overlap = fam_samples & anc_samples
            validation['sample_overlap'] = len(overlap)
            validation['fam_samples'] = len(fam_samples)

            if len(overlap) < len(fam_samples):
                missing = len(fam_samples) - len(overlap)
                validation['warnings'].append(f"{missing} genotype samples missing ancestry calls")

    except Exception as e:
        validation['warnings'].append(f"Could not check FAM file overlap: {str(e)}")

    # Final status
    if validation['errors']:
        validation['status'] = 'FAIL'
    elif validation['warnings']:
        validation['status'] = 'WARN'

    with open("ancestry_validation.json", 'w') as f:
        json.dump(validation, f, indent=2)

    print(f"Ancestry Validation: {validation['status']}")
    if 'ancestry_counts' in validation:
        print(f"  Total samples: {validation['n_samples']}")
        print("  Counts by ancestry:")
        for anc, cnt in validation['ancestry_counts'].items():
            print(f"    {anc}: {cnt}")

    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
    """
}

process CHECK_GWAS_TOOL_REQUIREMENTS {
    tag "$gwas_tool"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    val gwas_tool
    val has_kinship
    val has_local_ancestry
    val ancestry_group
    val is_binary_trait

    output:
    path "tool_requirements.json", emit: report
    env TOOL_OK, emit: tool_ok

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import json
    import os

    tool = "${gwas_tool}".lower()
    has_kinship = "${has_kinship}".lower() == 'true'
    has_la = "${has_local_ancestry}".lower() == 'true'
    ancestry = "${ancestry_group}"
    is_binary = "${is_binary_trait}".lower() == 'true'

    requirements = {
        'tool': tool,
        'ancestry': ancestry,
        'is_binary_trait': is_binary,
        'has_kinship': has_kinship,
        'has_local_ancestry': has_la,
        'status': 'OK',
        'warnings': [],
        'recommendations': []
    }

    admixed_groups = ['AAC', 'AHI', 'LAT1', 'LAT2', 'LATINO', 'ADMIXED']
    is_admixed = ancestry in admixed_groups

    # Tool-specific checks
    if tool == 'genesis':
        if not has_kinship:
            requirements['warnings'].append("GENESIS works best with a kinship matrix for related samples")
            requirements['recommendations'].append("Consider computing kinship matrix with PC-Relate or KING")

        if is_admixed:
            requirements['recommendations'].append("GENESIS is well-suited for admixed populations")

    elif tool == 'saige':
        if is_binary:
            requirements['recommendations'].append("SAIGE handles case-control imbalance well")
        if is_admixed:
            requirements['recommendations'].append("SAIGE is appropriate for admixed populations")

    elif tool == 'bolt-lmm':
        if is_admixed:
            requirements['warnings'].append("BOLT-LMM may not be optimal for highly admixed populations")
            requirements['recommendations'].append("Consider GENESIS or SAIGE for admixed samples")
        if is_binary:
            requirements['warnings'].append("BOLT-LMM is designed for quantitative traits")
            requirements['recommendations'].append("Consider SAIGE for binary traits")

    elif tool == 'regenie':
        requirements['recommendations'].append("REGENIE is computationally efficient for large datasets")

    elif tool == 'plink2':
        if has_kinship or is_admixed:
            requirements['warnings'].append("PLINK2 does not account for relatedness or population structure")
            requirements['recommendations'].append("Consider GENESIS or REGENIE for structured populations")

    # Tractor-specific checks
    if tool == 'tractor':
        if not has_la:
            requirements['status'] = 'FAIL'
            requirements['warnings'].append("Tractor requires local ancestry calls")
        if not is_admixed:
            requirements['warnings'].append("Tractor is designed for admixed populations")

    # Set environment variable for downstream use
    tool_ok = 'true' if requirements['status'] == 'OK' else 'false'

    with open("tool_requirements.json", 'w') as f:
        json.dump(requirements, f, indent=2)

    print(f"Tool: {tool}")
    print(f"Status: {requirements['status']}")
    if requirements['warnings']:
        print("Warnings:")
        for w in requirements['warnings']:
            print(f"  - {w}")
    if requirements['recommendations']:
        print("Recommendations:")
        for r in requirements['recommendations']:
            print(f"  - {r}")

    # Write env file
    with open('TOOL_OK', 'w') as f:
        f.write(tool_ok)
    """
}

process VALIDATE_FINE_MAPPING_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_reference

    output:
    tuple val(meta), path("${prefix}.fm_validation.json"), emit: report
    tuple val(meta), env(FM_READY), emit: ready

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import json
    import os
    import sys

    validation = {
        'status': 'READY',
        'warnings': [],
        'errors': [],
        'meta_id': '${meta.id}'
    }

    # Read summary statistics
    try:
        if "${sumstats}".endswith('.gz'):
            df = pd.read_csv("${sumstats}", sep='\\t', compression='gzip', nrows=1000)
        else:
            df = pd.read_csv("${sumstats}", sep='\\t', nrows=1000)

        # Check required columns for fine-mapping
        required = ['CHR', 'POS', 'SNP', 'BETA', 'SE']
        optional = ['P', 'MAF', 'N', 'Z']

        found_required = [c for c in required if c in df.columns]
        missing_required = [c for c in required if c not in df.columns]

        validation['found_columns'] = found_required
        validation['missing_required'] = missing_required

        if missing_required:
            validation['errors'].append(f"Missing required columns: {missing_required}")
            validation['status'] = 'NOT_READY'

        # Check for Z-score or compute from BETA/SE
        if 'Z' not in df.columns and 'BETA' in df.columns and 'SE' in df.columns:
            validation['warnings'].append("Z-scores will be computed from BETA/SE")

        # Check MAF availability
        if 'MAF' not in df.columns:
            validation['warnings'].append("MAF not found - some fine-mapping methods may not work optimally")

    except Exception as e:
        validation['errors'].append(f"Failed to read summary stats: {str(e)}")
        validation['status'] = 'NOT_READY'

    # Check LD reference
    ld_ref_path = "${ld_reference}"
    if not os.path.exists(ld_ref_path) or os.path.getsize(ld_ref_path) == 0:
        validation['warnings'].append("LD reference not provided or empty")
        validation['recommendations'] = ["Provide ancestry-matched LD reference for accurate fine-mapping"]

    # Save validation
    with open("${prefix}.fm_validation.json", 'w') as f:
        json.dump(validation, f, indent=2)

    fm_ready = 'true' if validation['status'] == 'READY' else 'false'
    print(f"Fine-mapping validation: {validation['status']}")

    # Write to env
    with open('FM_READY', 'w') as f:
        f.write(fm_ready)
    """
}

process DECIDE_FINE_MAPPING_METHOD {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(sumstats)
    val n_ancestries
    val n_traits
    val traits_related
    val default_within
    val default_multi

    output:
    tuple val(meta), env(FM_METHOD), emit: method

    script:
    """
    #!/usr/bin/env python3

    import os

    n_anc = int("${n_ancestries}")
    n_traits = int("${n_traits}")
    related = "${traits_related}".lower() == 'true'
    default_within = "${default_within}"
    default_multi = "${default_multi}"

    # Decision logic for fine-mapping method
    if n_anc == 1:
        # Single ancestry - use within-ancestry method
        method = default_within  # e.g., 'polyfun-susie'
    elif n_anc > 1:
        # Multi-ancestry
        if n_traits > 1 and related:
            # Multiple related traits - use MG-FLASH-FM
            method = 'mg-flash-fm'
        else:
            # Single trait or unrelated traits - use SuSIE-ME
            method = 'susie-me'
    else:
        method = default_within

    print(f"Selected fine-mapping method: {method}")
    print(f"  Ancestries: {n_anc}")
    print(f"  Traits: {n_traits}")
    print(f"  Traits related: {related}")

    with open('FM_METHOD', 'w') as f:
        f.write(method)
    """
}
