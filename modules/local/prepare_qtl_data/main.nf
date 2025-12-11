process PREPARE_QTL_DATA {
    tag "$name"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::pyyaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(name), val(type), path(qtl_file)

    output:
    tuple val(name), val(type), path("${name}.standardized.tsv.gz"), emit: standardized
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import gzip
    import sys

    # Read QTL data
    if "${qtl_file}".endswith('.gz'):
        df = pd.read_csv("${qtl_file}", sep='\\t', compression='gzip')
    else:
        df = pd.read_csv("${qtl_file}", sep='\\t')

    # Standardize column names
    col_mapping = {
        # Variant columns
        'variant_id': ['SNP', 'rsid', 'variant', 'snp_id', 'MarkerName'],
        'chr': ['CHR', 'chromosome', 'chrom', '#CHROM'],
        'pos': ['POS', 'position', 'BP', 'bp'],
        'ref': ['REF', 'ref', 'A2', 'other_allele'],
        'alt': ['ALT', 'alt', 'A1', 'effect_allele'],
        # Gene columns
        'gene_id': ['gene', 'gene_id', 'GENE', 'phenotype_id'],
        'gene_name': ['gene_name', 'gene_symbol', 'SYMBOL'],
        # Effect columns
        'beta': ['BETA', 'beta', 'effect', 'slope'],
        'se': ['SE', 'se', 'stderr', 'beta_se'],
        'pvalue': ['P', 'pvalue', 'p_value', 'pval', 'p.value'],
        'maf': ['MAF', 'maf', 'af', 'eaf'],
        'n': ['N', 'n_samples', 'sample_size']
    }

    standardized_cols = {}
    for new_col, old_cols in col_mapping.items():
        for old_col in old_cols:
            if old_col in df.columns:
                standardized_cols[old_col] = new_col
                break

    df = df.rename(columns=standardized_cols)

    # Add metadata columns
    df['qtl_type'] = '${type}'
    df['dataset'] = '${name}'

    # Ensure required columns exist
    required = ['variant_id', 'gene_id', 'beta', 'se', 'pvalue']
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Warning: Missing columns {missing}", file=sys.stderr)

    # Write standardized output
    df.to_csv("${name}.standardized.tsv.gz", sep='\\t', index=False, compression='gzip')

    print(f"Standardized {len(df)} QTL records for ${name}")

    # Write versions
    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "' + sys.version.split()[0] + '"\\n')
    """
}

process POOL_QTL_DATASETS {
    tag "$type"
    label 'process_medium'

    conda "conda-forge::python=3.11 conda-forge::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(type), val(names), path(qtl_files)

    output:
    tuple val(type), path("${type}.pooled.tsv.gz"), emit: pooled
    tuple val(type), path("${type}.pooled_stats.json"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import json
    import sys

    # Pool multiple QTL datasets of the same type
    # This creates a mega-set for colocalization across diverse populations

    qtl_files = "${qtl_files}".split()
    names = "${names}".split()

    pooled_dfs = []
    stats = {'type': '${type}', 'datasets': [], 'total_records': 0}

    for name, f in zip(names, qtl_files):
        try:
            df = pd.read_csv(f, sep='\\t', compression='gzip')
            df['source_dataset'] = name
            pooled_dfs.append(df)
            stats['datasets'].append({
                'name': name,
                'records': len(df),
                'genes': df['gene_id'].nunique() if 'gene_id' in df.columns else 0
            })
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)

    if pooled_dfs:
        pooled = pd.concat(pooled_dfs, ignore_index=True)

        # Remove duplicates (same variant-gene pair), keeping most significant
        if 'pvalue' in pooled.columns:
            pooled = pooled.sort_values('pvalue')
            pooled = pooled.drop_duplicates(
                subset=['variant_id', 'gene_id'],
                keep='first'
            )

        stats['total_records'] = len(pooled)
        stats['unique_genes'] = pooled['gene_id'].nunique() if 'gene_id' in pooled.columns else 0
        stats['unique_variants'] = pooled['variant_id'].nunique() if 'variant_id' in pooled.columns else 0

        pooled.to_csv("${type}.pooled.tsv.gz", sep='\\t', index=False, compression='gzip')
    else:
        # Create empty file if no data
        pd.DataFrame().to_csv("${type}.pooled.tsv.gz", sep='\\t', index=False, compression='gzip')

    with open("${type}.pooled_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)

    print(f"Pooled {stats['total_records']} QTL records from {len(stats['datasets'])} datasets")

    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "' + sys.version.split()[0] + '"\\n')
    """
}

process SELECT_QTL_BY_ANCESTRY {
    tag "$target_ancestry"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::pyyaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    val target_ancestry
    path qtl_config
    path available_datasets

    output:
    path "${target_ancestry}.selected_qtl_datasets.txt", emit: selected
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import yaml
    import sys

    # Load QTL configuration
    with open("${qtl_config}", 'r') as f:
        config = yaml.safe_load(f)

    # Get pooling strategy for target ancestry
    target = "${target_ancestry}"
    strategy = config.get('pooling_strategy', {}).get(target, {})

    primary = strategy.get('primary', [])
    secondary = strategy.get('secondary', [])

    # Load available datasets
    available = set()
    with open("${available_datasets}", 'r') as f:
        for line in f:
            available.add(line.strip())

    # Select datasets
    selected = []
    for ds in primary:
        if ds in available:
            selected.append(f"{ds}\\tprimary")
    for ds in secondary:
        if ds in available:
            selected.append(f"{ds}\\tsecondary")

    with open("${target_ancestry}.selected_qtl_datasets.txt", 'w') as f:
        f.write("\\n".join(selected))

    print(f"Selected {len(selected)} QTL datasets for {target}", file=sys.stderr)

    with open("versions.yml", 'w') as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "' + sys.version.split()[0] + '"\\n')
    """
}
