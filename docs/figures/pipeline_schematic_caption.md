# Pipeline schematic: figure caption and notes

Files in this directory:

| File | Use |
|------|-----|
| `pipeline_schematic.svg` / `.png` | Full-detail figure with Panels A, B and C. Sized for a 7.5 in page width; body text prints at about 6 pt. About 7.5 × 11 in, so use it as a full landscape page, scale to about 6.9 in wide, or drop Panel C. |
| `pipeline_schematic_AB_only.svg` / `.png` | Full-detail Panels A and B without the comparison table. About 7.5 × 9.5 in (fits one portrait page). |
| `pipeline_schematic_compact.svg` / `.png` | One line per module, Panels A, B and C. About 7.5 × 8.4 in. Recommended when the figure has to share a page with text. |
| `make_pipeline_schematic.py` (in the genotype pipeline repo, `documentation/figures/`) | Regenerates the SVGs (`python3 make_pipeline_schematic.py`, add `--compact` for the short version and `--no-table` to omit Panel C). No dependencies. Edit the `PANEL_A` / `PANEL_B` / `COMPARISON` lists to change wording. |
| `render_svg.js` (same location) | Renders an SVG to PNG with headless Chromium via Playwright: `node render_svg.js in.svg out.png 2` (the last number is the pixel scale; 2 gives about 400 dpi at 7.5 in). |

The SVG is editable in Illustrator, Inkscape, PowerPoint (Insert → Pictures) and can be uploaded to BioRender as an image asset.

---

## Figure caption (grant version)

**Figure 1. Modular Nextflow framework for ancestry-aware genotype harmonization and GWAS.**
**(A)** Genotype harmonization, imputation, quality control (QC) and ancestry inference (`Nextflow_Genotype_Pipeline`). Heterogeneous array genotypes (multiple platforms, batches, file layouts and genome builds) are converted, lifted over to hg38, aligned to the reference and merged by platform using a UNION strategy that preserves variant diversity (Module 1). Each platform is imputed separately on the TOPMed, Michigan and/or All of Us (AnVIL) imputation servers through automated job submission and retrieval (Module 2), filtered with the ancestry-aware MagicalRsq-X quality metric rather than a fixed R² cut-off (Module 3), harmonized and INTERSECTION-merged across platforms (Module 4), optionally re-imputed to recover variants lost at the intersection (Module 5), and passed through final sample- and variant-level QC including Hardy–Weinberg, minor allele frequency, heterozygosity, sex-check, relatedness and principal-component analyses (Module 6). Global ancestry is inferred with GRAF-anc and ADMIXTURE, and local ancestry with RFMix v2, FLARE or G-NOMIX against a 1000 Genomes–HGDP reference panel (Module 7). An optional benchmarking module compares imputed genotypes with whole-genome-sequencing truth, contrasts six QC/merge strategies across three imputation servers, and simulates phenotypes with GCTA to quantify GWAS power, hit recovery and local-ancestry accuracy within each ancestry group (Module 8; inset panels are illustrative). The output is an analysis-ready genotype set with global ancestry calls and local ancestry tracts.
**(B)** Ancestry-aware GWAS and downstream analyses (`Nextflow_GWAS_pipeline`). Harmonized genotypes, phenotypes and ancestry calls from (A) are QC'd and stratified into major and admixed ancestry groups, and GWAS is run in parallel per group and trait with REGENIE, SAIGE, BOLT-LMM, PLINK2 or GENESIS mixed models; admixed groups are additionally analysed with Tractor, which decomposes effects by local ancestry (two-way for African American, three-way for Latino participants), and time-to-event outcomes with SPA-Cox. Group-level results are combined with MR-MEGA meta-regression, which models ancestry-correlated heterogeneity, then fine-mapped within ancestry (PolyFun + SuSiE) and across ancestries (MG-FLASH-FM, SuSiE-ME). Downstream subworkflows compute ancestry-aware polygenic risk scores (PRS-CSx, PRS-CS, GAUDI, LDpred2), heritability and cross-ancestry genetic correlation (LDSC, GCTA), colocalization with ancestry-matched QTL resources, and functional annotation (MAGMA, FUMA, LAVA, FLAMES), followed by visualization and automated reporting. Every module in both pipelines is an independent Nextflow DSL2 entry point, so the framework can be run end-to-end or any module can be run alone on existing outputs; dashed outlines mark optional modules and orange marks admixed-population methods. **(C)** Principal departures from a typical GWAS workflow: end-to-end containerized merging instead of cohort-specific scripts; two-step imputation against the diverse TOPMed and All of Us AnVIL panels rather than a single largely European panel; ancestry-aware MagicalRsq-X filtering instead of a fixed R² cut-off; GRAF-anc plus local ancestry inference instead of self-report or global principal components; benchmarking built into the workflow; and containerized, parallelized, portable execution on HPC/SLURM.

---

## Shorter caption (about 120 words)

**Figure 1. Modular Nextflow framework for ancestry-aware genotype harmonization and GWAS.** **(A)** Array genotypes from multiple platforms and builds are harmonized, imputed on TOPMed, Michigan and/or All of Us servers, filtered with the ancestry-aware MagicalRsq-X metric, merged, QC'd and assigned global (GRAF-anc, ADMIXTURE) and local (RFMix, FLARE, G-NOMIX) ancestry. An optional module benchmarks imputation against WGS truth and simulates GWAS power by ancestry (insets illustrative). **(B)** Harmonized genotypes and ancestry calls feed ancestry-stratified GWAS (REGENIE, SAIGE, GENESIS; Tractor for admixed groups; SPA-Cox for survival), MR-MEGA meta-analysis, within- and multi-ancestry fine-mapping, and downstream PRS, heritability, colocalization and functional annotation with reporting. **(C)** Departures from a typical GWAS workflow: containerized end-to-end merging, two-step imputation on diverse panels, ancestry-aware variant filtering, GRAF-anc plus local ancestry, built-in benchmarking and portable execution. Modules run end-to-end or individually; dashed outlines are optional.

---

## Notes for adapting the figure

- The five inset charts in Module 8 are schematic placeholders that show what the benchmarking module reports (concordance with WGS truth, variant retention under MagicalRsq-X versus a fixed R² filter, imputation R² by ancestry, simulated GWAS power for approaches A–F, and local-ancestry accuracy). Replace them with real panels from `Module8_Benchmarking.nf` output once benchmark runs are complete, and drop the "illustrative" footnote.
- Approach letters A–F refer to the six QC/merge strategies in `documentation/BENCHMARKING_APPROACH_WORKFLOWS.md` (A–D: traditional R² filtering with different QC/merge timing; E–F: this pipeline with one- or two-step MagicalRsq-X filtering).
- Ancestry group codes follow GRAF-anc: EUR, AFR, EAS, SAS, AMR, MID, plus admixed AAC (African American), AHI/LAT1/LAT2 (Latino).
