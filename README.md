# Integrated Nipah Virus Host-Response Transcriptomics

**GitHub repository:** https://github.com/iamsuwaib/Integrated-Nipah-Virus-Host-Response-Transcriptomics

**Associated manuscript:**  
Khan S, Khan TU, Ahmad J, Ullah A, Butt S, Ahmed A. *Integrated host-response transcriptomics separates a conserved antiviral-associated program from in vivo tissue-associated disease programs in Nipah virus infection.* Computational and Structural Biotechnology Reports (under review), 2026.

---

## Overview

This repository contains the R and Python analysis scripts used to generate all figures, tables, and supplementary outputs for the manuscript. The analysis integrates three public GEO series to define conserved antiviral and tissue-level disease-amplification programs during Nipah virus (NiV) infection:

| Dataset | System | Platform | Comparison |
|---|---|---|---|
| GSE32902 | Human HUVEC | Microarray | NiV vs. mock |
| GSE33133 | Human HUVEC | Microarray | NiV, NiV-dC vs. mock (the NiV vs. mock contrast re-lists the GSE32902 sample records) |
| GSE310471 | AGM lung & tonsil | RNA-seq | Baseline, 3, 4, 5 DPI |

GSE32902 and GSE33133 are two related GEO series representing one underlying HUVEC experiment (the GSE33133 NiV vs. mock contrast re-uses the four GSE32902 sample records and is excluded from recurrence counting and axis means), whereas GSE310471 is an independent in vivo AGM study. Study-level and leave-one-system-out sensitivity analyses (scripts 16 and 21) treat these as two independent evidence units.

Data are available from the NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/).

---

## Analysis Scripts

Scripts are numbered in the order they should be run (script 25 requires the network objects from script 07, and scripts 26 and 27 are run together to produce the orthology table). Each script is self-contained and outputs figures and/or result tables.

| Script | Analysis |
|---|---|
| `01_GSE32902_limma_analysis.R` | Differential expression of GSE32902 HUVEC microarray (limma) |
| `02_GSE33133_limma_analysis_png_only.R` | Differential expression of GSE33133 HUVEC microarray, WT NiV and NiV-dC (limma) |
| `03_GSE310471_DESeq2_analysis_png_only.R` | Differential expression of GSE310471 AGM lung and tonsil RNA-seq (DESeq2); HUVEC-signature projection into tissues (Figure 3) |
| `04_integrated_cross_dataset_signature_heatmap.R` | Cross-dataset signature integration and recurrence counting (GSE33133 NiV vs. mock displayed with an asterisk and excluded from scoring); writes `integrated_signature_long.csv` (Table S4) |
| `05_candidate_biomarker_table.R` | Candidate shortlist and Table 1 (candidate host-response signatures; priority score defined from total, HUVEC-specific and in vivo-specific recurrence) |
| `06_tf_upstream_regulator_signature_scoring.R` | Targeted regulatory-signature scoring for six axes (Figure 4B, Table S6); scoring function shared with script 20 through `regulator_axis_common.R` |
| `07_wgcna_gse310471_exploratory.R` | WGCNA co-expression network analysis for GSE310471 lung and tonsil (Figures 4A, S16, S17) |
| `08_secretome_surfaceome.R` | Extracellular and cell-surface mediator prioritization (Figure 6, Tables S13-S15); evidence tiers count GSE32902 and GSE33133 as one HUVEC evidence unit |
| `09_wgcna_candidate_signature_membership.R` | Candidate gene WGCNA module membership mapping (Figure S18, Tables S7-S8) |
| `10_focused_wgcna_module_discovery_figure.R` | Focused WGCNA module discovery figure (Figure 4A) |
| `11_gsea_module_enrichment.R` | Ranked GO Biological Process GSEA and module ORA (Figure 5, Tables S11-S12) |
| `12_targeted_gsea_figure_refinement.R` | Targeted GSEA and module enrichment figure refinement (Figure 5) |
| `13_combine_figure3_panels.py` | Figure 3 panel compositing (lung/tonsil projection heatmaps) |
| `14_combine_figure2_panels.py` | Earlier Figure 2 panel compositing (retained for reference; superseded by script 28 for the print-scale Figure 2A/2B) |
| `15_combine_figure4_panels.py` | Earlier Figure 4 panel compositing (retained for reference; superseded by script 29) |
| `16_sensitivity_rank_based_scoring.R` | Sensitivity analysis: priority score vs. an alternative rank-based composite score (Supplementary Figure S20, Tables S16-S17) |
| `17_cell_composition_marker_scores.R` | Cell-composition marker-gene scoring across infection time using gene-wise z-scored marker sets, Benjamini-Hochberg correction across marker-set correlations, and a broader lung-specific marker panel in tonsil (Supplementary Figures S21-S22, Tables S18-S19) |
| `18_progression_lrt_and_trajectories.R` | DESeq2 likelihood-ratio test for an infection-time effect (DPI as a categorical factor); module eigengene trajectories (Table S20, Supplementary Figure S23) |
| `19_wgcna_module_stability.R` | Bootstrap and leave-one-out WGCNA module-trait correlation stability, and full network re-clustering (Table S21, Supplementary Figure S26) |
| `20_regulator_axis_permutation_test.R` | Genome-wide permutation-null testing of the regulatory-signature axis scores: one common random gene set per replicate evaluated across all contrasts (primary design), with an independent per-contrast draw as a sensitivity design; 10,000 replicates; Benjamini-Hochberg and Holm adjustment across the six axes (Table S22, Supplementary Figure S27) |
| `21_study_level_sensitivity.R` | Equal-weight-by-study composite score and leave-one-system-out stability of the candidate shortlist (Tables S24-S25, Supplementary Figures S28-S29) |
| `22_tonsil_module_alveolar_gene_robustness.R` | Tonsil blue/ME2 module recomputed after removing lung-specific genes (Table S26, Supplementary Figure S30) |
| `23_tonsil_lung_admixture_model.R` | Gene-level lung-RNA admixture model for the tonsil infection response (Table S27, Supplementary Figure S31) |
| `24_ordered_time_trend_tests.R` | Ordered-time (numeric DPI) trend tests for module eigengenes and genes (Table S28, Supplementary Figure S32) |
| `25_wgcna_soft_power_sensitivity.R` | Explicit soft-power rule (pickSoftThreshold estimate at R2 >= 0.85, otherwise fixed fallback of 6) and sensitivity of the lung tan/ME12 and tonsil blue/ME2 modules across nearby powers (Table S29, Supplementary Figures S24-S25 and S33) |
| `26_orthology_verification.R` | Verification of human to African green monkey one-to-one orthology for the final candidate set using Ensembl Compara (requires internet access) (Table S30) |
| `27_orthology_verification_ncbi.R` | Orthology verification with NCBI Gene orthologs and E-utilities, used to resolve genes not returned by Ensembl; combined with script 26 output (Table S30) |
| `28_print_scale_figures.py` | Python (matplotlib) rebuild of Figure 2A (integrated heatmap) and Figure 2B (conserved-signature lollipop) at print size (6.5 in wide, 600 dpi) so that labels remain legible at page width |
| `29_print_scale_figure4.py` | Python (matplotlib) rebuild of Figure 4A (WGCNA module membership panels) and Figure 4B (regulatory-signature heatmap) at print size |
| `31_print_scale_figure5.py` | Python (matplotlib) rebuild of Figure 5A (targeted GSEA dot plot) and Figure 5B (module over-representation) at print size from Supplementary Tables S11 and S12 |
| `32_print_scale_figure6.py` | Python (matplotlib) rebuild of Figure 6A (candidate ranking by study-level priority score) and Figure 6B (mediator classes) at print size from the refined secretome/surfaceome tables |
| `30_workflow_schematic.py` | Python (matplotlib) rebuild of Figure 1 (study design and workflow schematic) |
| `regulator_axis_common.R` | Shared helper defining the six regulatory axes, the scoring function (mean log2FC of detected axis genes x log2(n detected + 1)) and the permutation universe; sourced by scripts 06 and 20 so that both use one authoritative implementation |

---

## Requirements

R version 4.4.2. Key packages:

- `limma` v3.62.2
- `DESeq2` v1.46.0
- `WGCNA` (for scripts 07, 09, 10, 19, 22, 24, 25)
- `clusterProfiler` v4.14.6
- `fgsea` v1.32.4
- `ggplot2` v3.5.1
- `pheatmap` v1.0.13
- `dplyr` v1.1.4
- `tidyverse` v2.0.0

Python scripts (13 to 15, 28 to 32) require Python 3.9 or later with `numpy`, `pandas`, `scipy`, `matplotlib` and `Pillow`. Scripts 26 and 27 query the Ensembl REST service and NCBI E-utilities and therefore require internet access.

Install Bioconductor packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "DESeq2", "clusterProfiler", "fgsea", "WGCNA"))
```

---

## Data Availability

All raw data are publicly available from NCBI GEO. No patient data or restricted datasets are included in this repository. Processed supplementary tables (S1-S30) and supplementary figures (S1-S33) are provided with the manuscript submission.

---

## Citation

If you use these scripts, please cite the associated manuscript (full citation to be updated upon acceptance) and, for reproducibility, the specific tagged release of this repository used (see Releases).

## Contact

Sohaib Khan — sohaib.khan@ug.edu.pl  
International Centre for Cancer Vaccine Science, University of Gdansk, Poland
