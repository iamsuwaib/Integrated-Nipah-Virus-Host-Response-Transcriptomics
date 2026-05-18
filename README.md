# Nipah Virus Host-Response Transcriptomics

This repository contains the R analysis scripts used for the manuscript:

**A Two-Layer Transcriptomic Model of Nipah Virus Infection Identifies Conserved Antiviral Defense and In Vivo Tissue Disease Programs**

The analysis integrates public GEO transcriptomic datasets from human endothelial-cell Nipah virus infection and African green monkey lung and tonsil infection to identify conserved antiviral/interferon responses and in vivo tissue-associated complement, coagulation/hemostasis, vascular-inflammatory, and metabolic-remodeling programs.

## Public Datasets

| Dataset | System | Platform type | Main comparison |
|---|---|---|---|
| GSE32902 | Human HUVEC | Expression array | NiV vs mock |
| GSE33133 | Human HUVEC | Expression array | NiV, NiV-dC, mock |
| GSE310471 | African green monkey lung and tonsil | RNA-seq counts | 3, 4, and 5 DPI vs baseline |

GSE33133 extends the HUVEC endothelial infection design by adding the NiV-dC condition; it should not be treated as a fully independent validation dataset for GSE32902. GSE310471 provides the independent in vivo progression layer.

## Analysis Overview

The workflow includes:

1. Dataset-level quality control and differential expression analysis.
2. Cross-dataset integration at the gene-symbol and contrast-result level.
3. Candidate biomarker prioritization.
4. Targeted regulator-axis scoring.
5. Exploratory WGCNA on GSE310471 lung and tonsil expression matrices.
6. Ranked GO Biological Process GSEA and WGCNA module-enrichment analysis.
7. Targeted manuscript figure generation from reviewed analysis outputs.

Raw expression values were not pooled across platforms, species, or biological systems. Each dataset was normalized and analyzed independently using dataset-appropriate statistical methods, and integration was performed using contrast-level statistics.

## Script Order

Run scripts from the project root after downloading and extracting the required GEO files into the expected local folders.

| Script | Purpose |
|---|---|
| `01_GSE32902_limma_analysis.R` | QC and limma differential expression for GSE32902 |
| `02_GSE33133_limma_analysis_png_only.R` | QC and limma differential expression for GSE33133 |
| `03_GSE310471_DESeq2_analysis_png_only.R` | QC and DESeq2 differential expression for GSE310471 |
| `04_integrated_cross_dataset_signature_heatmap.R` | Integrated cross-dataset signature matrix and heatmap |
| `05_candidate_biomarker_table.R` | Candidate biomarker/signature prioritization |
| `06_tf_upstream_regulator_signature_scoring.R` | Targeted regulator-axis scoring |
| `07_wgcna_gse310471_exploratory.R` | Exploratory WGCNA for lung and tonsil |
| `08_extract_wgcna_hubs_from_existing_rds.R` | Extract module hub/candidate information from WGCNA objects |
| `09_wgcna_candidate_signature_membership.R` | Candidate gene module-membership analysis |
| `10_focused_wgcna_module_discovery_figure.R` | Focused WGCNA manuscript figure generation |
| `14_gsea_module_enrichment.R` | Ranked GO BP GSEA and WGCNA module ORA |
| `15_targeted_gsea_figure_refinement.R` | Targeted GSEA/module-enrichment figure refinement |

Manuscript assembly, exploratory drug ranking, and unused validation-triage helper scripts are intentionally excluded from this final GitHub-ready folder.

## Main R Dependencies

The manuscript analyses were performed in R v4.4.2. Major packages included:

- `limma` v3.62.2
- `DESeq2` v1.46.0
- `clusterProfiler` v4.14.6
- `fgsea` v1.32.4
- `WGCNA`
- `ggplot2` v4.0.2
- `pheatmap` v1.0.13
- `dplyr` v1.1.4
- `tidyverse` v2.0.0

Package versions may vary across systems; scripts include package-loading sections and should be checked against the accompanying manuscript methods.

## Outputs

Scripts generate dataset-level results, integrated tables, candidate biomarker tables, WGCNA module outputs, GSEA/module-enrichment tables, and PNG figures. Output folders are defined inside each script and may need adjustment if the project is moved to a different directory.

## Citation

If using these scripts, please cite the associated manuscript and the original GEO datasets.

