# Integrated Nipah Virus Host-Response Transcriptomics

**GitHub repository:** https://github.com/iamsuwaib/Integrated-Nipah-Virus-Host-Response-Transcriptomics

**Associated manuscript:**  
Khan S, Khan TU, Ahmad J, Ullah A, Butt S, Ahmed A. *Integrated host-response transcriptomics separates a conserved antiviral-associated program from in vivo tissue-associated disease programs in Nipah virus infection.* Computational and Structural Biotechnology Reports (under revision), 2026.

---

## Overview

This repository contains the R analysis scripts used to generate all figures, tables, and supplementary outputs for the manuscript. The analysis integrates three public GEO datasets to define conserved antiviral and tissue-level disease-amplification programs during Nipah virus (NiV) infection:

| Dataset | System | Platform | Comparison |
|---|---|---|---|
| GSE32902 | Human HUVEC | Microarray | NiV vs. mock |
| GSE33133 | Human HUVEC | Microarray | NiV, NiV-dC vs. mock |
| GSE310471 | AGM lung & tonsil | RNA-seq | Baseline, 3, 4, 5 DPI |

Data are available from the NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/).

---

## Analysis Scripts

Scripts are numbered in the order they should be run. Each script is self-contained and outputs figures and/or result tables.

| Script | Analysis |
|---|---|
| `01_GSE32902_limma_analysis.R` | Differential expression of GSE32902 HUVEC microarray (limma) |
| `02_GSE33133_limma_analysis_png_only.R` | Differential expression of GSE33133 HUVEC microarray, WT NiV and NiV-dC (limma) |
| `03_GSE310471_DESeq2_analysis_png_only.R` | Differential expression of GSE310471 AGM lung and tonsil RNA-seq (DESeq2); HUVEC-signature projection into tissues (Figure 3) |
| `04_integrated_cross_dataset_signature_heatmap.R` | Cross-dataset signature integration and heatmap (Figure 2A, 2B) |
| `05_candidate_biomarker_table.R` | Candidate biomarker shortlist and Table 1 |
| `06_tf_upstream_regulator_signature_scoring.R` | Targeted regulatory-signature scoring (Figure 4B, Table S6) |
| `07_wgcna_gse310471_exploratory.R` | WGCNA co-expression network analysis for GSE310471 lung and tonsil (Figures 4A, S16, S17) |
| `08_secretome_surfaceome.R` | Extracellular and cell-surface mediator prioritization (Figure 6, Tables S13-S15) |
| `09_wgcna_candidate_signature_membership.R` | Candidate gene WGCNA module membership mapping (Figure S18, Tables S7-S8) |
| `10_focused_wgcna_module_discovery_figure.R` | Focused WGCNA module discovery figure (Figure 4A) |
| `11_gsea_module_enrichment.R` | Ranked GO Biological Process GSEA and module ORA (Figure 5, Tables S11-S12) |
| `12_targeted_gsea_figure_refinement.R` | Targeted GSEA and module enrichment figure refinement (Figure 5) |
| `13_combine_figure3_panels.py` | Figure 3 panel compositing (lung/tonsil projection heatmaps) |
| `14_combine_figure2_panels.py` | Figure 2 panel compositing (integrated heatmap + conserved-signature lollipop) |
| `15_combine_figure4_panels.py` | Figure 4 panel compositing (WGCNA module panel + regulatory-signature heatmap) |
| `16_sensitivity_rank_based_scoring.R` | Sensitivity analysis: priority score vs. an independent rank-based composite score (Supplementary Figure S20, Tables S16-S17) |
| `17_cell_composition_marker_scores.R` | Cell-composition marker-gene scoring across infection time (Supplementary Figures S21-S22, Tables S18-S19) |
| `18_progression_lrt_and_trajectories.R` | Formal DESeq2 likelihood-ratio test for the infection-time effect; module eigengene trajectories (Table S20, Supplementary Figure S23) |
| `19_wgcna_module_stability.R` | Bootstrap and leave-one-out WGCNA module-trait correlation stability, and full network re-clustering (Table S21, Supplementary Figure S26) |
| `20_regulator_axis_permutation_test.R` | Genome-wide permutation-null testing of the regulatory-signature axis scores (Table S22, Supplementary Figure S27) |

---

## Requirements

R version 4.4.2. Key packages:

- `limma` v3.62.2
- `DESeq2` v1.46.0
- `WGCNA` (for scripts 07, 09, 10, 19)
- `clusterProfiler` v4.14.6
- `fgsea` v1.32.4
- `ggplot2` v3.5.1
- `pheatmap` v1.0.13
- `dplyr` v1.1.4
- `tidyverse` v2.0.0

Install Bioconductor packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "DESeq2", "clusterProfiler", "fgsea", "WGCNA"))
```

---

## Data Availability

All raw data are publicly available from NCBI GEO. No patient data or restricted datasets are included in this repository. Processed supplementary tables (S1-S23) and supplementary figures (S1-S27) are provided with the manuscript submission.

---

## Citation

If you use these scripts, please cite the associated manuscript (full citation to be updated upon acceptance) and, for reproducibility, the specific tagged release of this repository used (see Releases).

## Contact

Sohaib Khan — sohaib.khan@ug.edu.pl  
International Centre for Cancer Vaccine Science, University of Gdansk, Poland
