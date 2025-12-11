# Ancestry-Aware GWAS Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A comprehensive Nextflow pipeline for multi-ancestry genome-wide association studies (GWAS) with fine-mapping, meta-analysis, polygenic risk score calculation, and functional annotation.

## Overview

This pipeline addresses the critical need for ancestry-aware genetic analysis by:

1. **Ancestry-stratified GWAS**: Runs separate GWAS for each major ancestry group (EUR, AFR, EAS, SAS, AMR, MID, AAC, AHI, HET)
2. **Multi-ancestry meta-analysis**: Uses MR-MEGA to model ancestry-correlated heterogeneity
3. **Within-ancestry fine-mapping**: PolyFun+SuSIE for annotation-prioritized fine-mapping
4. **Multi-ancestry fine-mapping**: MG-FLASH-FM for related traits, SuSIE-ME for single traits
5. **Ancestry-aware PRS**: Multiple methods including PRS-CSx, GAUDI, and local ancestry-informed approaches
6. **Heritability analysis**: Ancestry-specific h² and cross-ancestry genetic correlations
7. **Functional annotation**: MAGMA, FUMA, LAVA, FLAMES with cell-type specificity

## Pipeline Architecture

```
┌─────────────────┐
│   Input Data    │
│  (Genotypes +   │
│   Phenotypes)   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   QC + Ancestry │
│    Inference    │
│   (GRAF-ANC)    │
└────────┬────────┘
         │
    ┌────┴────┬────────┬────────┬────────┬────────┬────────┬────────┐
    ▼         ▼        ▼        ▼        ▼        ▼        ▼        ▼
┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐
│  EUR  │ │  AFR  │ │  EAS  │ │  SAS  │ │  AMR  │ │  MID  │ │  AAC  │ │  AHI  │
│ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │
└───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘
    │         │        │        │        │        │        │        │
    └────┬────┴────────┴────────┴────────┴────────┴────────┴────────┘
         │
         ▼
┌─────────────────┐
│   MR-MEGA       │
│  Meta-Analysis  │
└────────┬────────┘
         │
    ┌────┴─────────────┐
    ▼                  ▼
┌─────────────┐  ┌─────────────┐
│ Within-Anc  │  │ Multi-Anc   │
│ Fine-Map    │  │ Fine-Map    │
│ PolyFun+    │  │ MG-FLASH-FM │
│ SuSIE       │  │ or SuSIE-ME │
└──────┬──────┘  └──────┬──────┘
       │                │
       └────────┬───────┘
                ▼
┌───────────────────────────────────────────────────────────┐
│                    Downstream Analysis                     │
├─────────────┬─────────────┬──────────────┬────────────────┤
│     PRS     │ Heritability│ Colocalization│   Functional   │
│  PRS-CSx    │    LDSC     │    coloc     │     MAGMA      │
│  PRS-CS     │    GCTA     │   eCAVIAR    │     FUMA       │
│   GAUDI     │  Cross-Anc  │  fastENLOC   │     LAVA       │
│  LDpred2    │     rg      │              │    FLAMES      │
└─────────────┴─────────────┴──────────────┴────────────────┘
```

## Quick Start

```bash
# Basic run with REGENIE
nextflow run main.nf \
  --input samplesheet.csv \
  --phenotype_file phenotypes.txt \
  --phenotype_cols "trait1,trait2" \
  --outdir results \
  -profile docker

# With all analyses enabled
nextflow run main.nf \
  --input samplesheet.csv \
  --phenotype_file phenotypes.txt \
  --phenotype_cols "BMI,Height,Weight" \
  --covariate_cols "age,sex,PC1,PC2,PC3,PC4,PC5" \
  --gwas_tool regenie \
  --meta_analysis true \
  --fine_mapping true \
  --prs_analysis true \
  --prs_methods "prs-csx,gaudi,ldpred2" \
  --heritability true \
  --functional_annotation true \
  --outdir results \
  -profile singularity
```

## Input Requirements

### Samplesheet

Create a CSV file with genotype file paths:

```csv
sample_id,cohort,bed,bim,fam
study1,cohort_A,/path/to/study1.bed,/path/to/study1.bim,/path/to/study1.fam
study2,cohort_B,/path/to/study2.bed,/path/to/study2.bim,/path/to/study2.fam
```

Supported formats: PLINK1 (bed/bim/fam), PLINK2 (pgen/pvar/psam), BGEN, VCF

### Phenotype File

Tab-separated file with sample IDs and phenotype values:

```
IID	trait1	trait2	age	sex	PC1	PC2	PC3	PC4	PC5
sample1	25.3	170.2	45	1	0.01	-0.02	0.005	0.001	-0.003
sample2	28.1	165.8	38	2	-0.01	0.03	-0.002	0.004	0.001
```

## Key Parameters

### Ancestry Inference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ancestry_method` | `graf` | Method: graf, admixture, pca |
| `--ancestry_groups` | `EUR,AFR,EAS,SAS,AMR,MID,AAC,AHI,HET` | Ancestry groups to analyze |
| `--min_ancestry_n` | `100` | Minimum samples per ancestry |
| `--local_ancestry` | `false` | Run local ancestry inference |

### GWAS

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gwas_tool` | `regenie` | Tool: regenie, saige, bolt-lmm, plink2 |
| `--gwas_model` | `additive` | Model: additive, dominant, recessive |
| `--p_threshold` | `5e-8` | Genome-wide significance threshold |

### Meta-Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--meta_method` | `mr-mega` | Method: mr-mega, metal, metasoft |
| `--mr_mega_env` | `false` | Include environmental variable |

### Fine-Mapping

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--within_ancestry_fm` | `polyfun-susie` | Within-ancestry method |
| `--multi_ancestry_fm` | `mg-flash-fm` | Multi-ancestry method for related traits |
| `--single_trait_fm` | `susie-me` | Multi-ancestry method for single traits |
| `--related_traits` | `true` | Whether traits are related |

### PRS

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--prs_methods` | `prs-csx,prs-cs,gaudi,ldpred2` | PRS methods to run |
| `--prs_validation` | `true` | Validate PRS in target cohort |
| `--prs_best_method` | `true` | Determine best method per ancestry |

## Outputs

```
results/
├── ancestry/
│   ├── graf/                    # GRAF ancestry calls
│   └── stratified/              # Ancestry-stratified genotypes
├── gwas/
│   └── {ancestry}/{trait}/      # Per-ancestry GWAS results
├── meta_analysis/
│   └── {trait}/                 # MR-MEGA meta-analysis
├── fine_mapping/
│   ├── {ancestry}/{trait}/      # Within-ancestry fine-mapping
│   └── multi_ancestry/{trait}/  # Multi-ancestry fine-mapping
├── prs/
│   └── {trait}/                 # PRS weights and scores
├── heritability/
│   └── {ancestry}/{trait}/      # h² estimates
├── functional/
│   └── {trait}/                 # Functional annotations
├── colocalization/
│   └── {trait}/                 # Colocalization results
├── plots/
│   └── {trait}/                 # Manhattan, QQ, regional plots
└── reports/
    ├── multiqc_report.html      # MultiQC report
    └── ancestry_gwas_report.html # Custom pipeline report
```

## Ancestry Groups (GRAF-ANC)

| Code | Description |
|------|-------------|
| EUR | European |
| AFR | African |
| EAS | East Asian |
| SAS | South Asian |
| AMR | Native American |
| MID | Middle Eastern |
| AAC | African American (admixed) |
| AHI | American Hispanic (admixed) |
| HET | Heterogeneous/Other |

## Methods Description

### MR-MEGA Meta-Analysis
MR-MEGA performs meta-regression of genetic association data, modeling allelic effects as a function of axes of genetic variation. This accounts for:
- Ancestry-correlated heterogeneity (differences due to genetic ancestry)
- Residual heterogeneity (other sources of between-study variation)

### PolyFun+SuSIE Fine-Mapping
Combines functional annotation-informed prior probabilities (PolyFun) with Sum of Single Effects regression (SuSIE) for:
- Annotation-prioritized causal variant identification
- Multiple causal variant detection per locus
- 95% credible set construction

### MG-FLASH-FM
Multi-trait multi-ancestry fine-mapping designed for related traits:
- Leverages shared genetic architecture across traits
- Accounts for ancestry-specific LD patterns
- Improves power through trait correlation

### PRS-CSx
Cross-population polygenic prediction using coupled continuous shrinkage priors:
- Jointly models summary statistics from multiple populations
- Shares information across ancestries for improved accuracy
- Particularly effective for underrepresented populations

### GAUDI
Local ancestry-informed PRS for admixed individuals:
- Uses local ancestry tracts to weight PRS contributions
- Addresses ancestry-specific effect sizes
- Optimal for admixed populations like African Americans

## Citation

If you use this pipeline, please cite:

```
[Pipeline citation - to be added]

Key methods:
- REGENIE: Mbatchou et al. (2021) Nature Genetics
- MR-MEGA: Magi et al. (2017) Nature Genetics
- PolyFun: Weissbrod et al. (2020) Nature Genetics
- SuSIE: Wang et al. (2020) J R Stat Soc
- MG-FLASH-FM: [Citation]
- PRS-CSx: Ruan et al. (2022) Nature Genetics
- GAUDI: [Citation]
```

## Requirements

- Nextflow >= 23.04.0
- Container runtime: Docker, Singularity, or Apptainer
- OR Conda/Mamba for local execution

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please read our [contributing guidelines](CONTRIBUTING.md) first.

## Support

- [Open an issue](https://github.com/your-org/ancestry-aware-gwas/issues)
- [Documentation](https://your-org.github.io/ancestry-aware-gwas/)
