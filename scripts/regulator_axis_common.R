###############################################################################
# Shared definitions for the targeted regulatory-signature axis analyses
# (used by 06_tf_upstream_regulator_signature_scoring.R and
# 20_regulator_axis_permutation_test.R)
#
# Why this file exists:
#   Scripts 06 and 20 share one set of axis definitions and read the same gene
#   universe, the genome-wide per-contrast differential-expression tables, so
#   every axis has one authoritative score. Axis genes need not be members of
#   the 57-gene curated integrated-signature panel (integrated_signature_long.csv);
#   for example the NF-kB / inflammatory chemokine targets NFKB1, RELA and
#   NFKBIA are not in the panel but are scored whenever they are detected in a
#   contrast.
#
# Scoring formula: score = mean(log2FC of detected axis genes) *
#   log2(number of detected axis genes + 1), per contrast; the axis mean score
#   is the average over the 9 scored contrasts (GSE33133 NiV-vs-mock is
#   displayed but excluded from scoring because it duplicates the GSE32902
#   sample records).
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"

tf_sets <- tribble(
  ~regulator, ~target_genes,
  "IRF7 / antiviral IRF axis", "MX1,MX2,OAS1,OAS2,OAS3,OASL,IFIT1,IFIT2,IFIT3,IFIT5,IFIH1,DDX58,RSAD2,ISG15,USP18,HERC5,HERC6,PARP9",
  "ISGF3-like STAT1/STAT2/IRF9 axis", "STAT1,STAT2,IRF9,MX1,MX2,OAS1,OAS2,OAS3,OASL,IFIT1,IFIT2,IFIT3,ISG15,USP18,RSAD2",
  "RIG-I/MDA5 sensing axis", "DDX58,IFIH1,IRF7,IRF9,STAT1,STAT2,IFIT1,IFIT2,IFIT3,MX1,OAS1,CXCL10,CXCL11",
  "NF-kB / inflammatory chemokine axis", "NFKB1,RELA,NFKBIA,TNF,IL6,CXCL10,CXCL11,CXCL9,CCL2,CCL5,ICAM1,VCAM1,SELE",
  "Endothelial activation axis", "ICAM1,VCAM1,SELE,ANGPT2,VWF,THBD,SERPINE1,PLAU,PLAUR,F3",
  "Complement/coagulation axis", "C3,C4A,C4B,C1QA,C1QB,C1QC,CFB,CFD,CFH,CFI,SERPING1,FGB,FGG,FGA,PROC,PROS1,THBD,SERPINE1"
) %>%
  mutate(target_genes = str_split(target_genes, ",")) %>%
  tidyr::unnest(target_genes) %>%
  mutate(target_genes = str_trim(target_genes))

axes <- unique(tf_sets$regulator)

# Differential-expression thresholds used throughout the study
SIG_PADJ <- 0.05
SIG_LOG2FC <- 1

read_limma_universe <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(SYMBOL = na_if(str_trim(as.character(SYMBOL)), "")) %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(SYMBOL, adj.P.Val, desc(abs(logFC))) %>%
    group_by(SYMBOL) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene = SYMBOL, log2FC = as.numeric(logFC), padj = as.numeric(adj.P.Val))
}

read_deseq2_universe <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(log2FoldChange)) %>%
    transmute(gene = as.character(Gene), log2FC = as.numeric(log2FoldChange), padj = as.numeric(padj)) %>%
    distinct(gene, .keep_all = TRUE)
}

contrast_registry <- tribble(
  ~contrast, ~file, ~type, ~scoring_included,
  "GSE32902_HUVEC_NiV_vs_Mock", "GSE32902/results/tables/GSE32902_limma_all_probes_annotated.csv", "limma", TRUE,
  "GSE33133_HUVEC_NiV_vs_Mock", "GSE33133/results/tables/GSE33133_limma_all_NiV_vs_Mock_annotated.csv", "limma", FALSE,
  "GSE33133_HUVEC_NiVdC_vs_Mock", "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_Mock_annotated.csv", "limma", TRUE,
  "GSE33133_HUVEC_NiVdC_vs_NiV", "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_NiV_annotated.csv", "limma", TRUE,
  "GSE310471_Lung_3DPI_vs_baseline", "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_3DPI_vs_baseline.csv", "deseq2", TRUE,
  "GSE310471_Lung_4DPI_vs_baseline", "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_4DPI_vs_baseline.csv", "deseq2", TRUE,
  "GSE310471_Lung_5DPI_vs_baseline", "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_5DPI_vs_baseline.csv", "deseq2", TRUE,
  "GSE310471_Tonsil_3DPI_vs_baseline", "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_3DPI_vs_baseline.csv", "deseq2", TRUE,
  "GSE310471_Tonsil_4DPI_vs_baseline", "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_4DPI_vs_baseline.csv", "deseq2", TRUE,
  "GSE310471_Tonsil_5DPI_vs_baseline", "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_5DPI_vs_baseline.csv", "deseq2", TRUE
)

load_universes <- function(scored_only = TRUE) {
  reg <- if (scored_only) dplyr::filter(contrast_registry, scoring_included) else contrast_registry
  out <- vector("list", nrow(reg))
  names(out) <- reg$contrast
  for (i in seq_len(nrow(reg))) {
    path <- file.path(project_dir, reg$file[i])
    if (!file.exists(path)) stop("Missing genome-wide table: ", path)
    out[[i]] <- if (reg$type[i] == "limma") read_limma_universe(path) else read_deseq2_universe(path)
  }
  out
}

# Score one gene set in one contrast universe (gene, log2FC, padj)
score_geneset <- function(universe_df, genes) {
  sub <- universe_df[universe_df$gene %in% genes, , drop = FALSE]
  n <- nrow(sub)
  if (n == 0) {
    return(tibble(n_targets_detected = 0L, n_targets_sig_up = NA_integer_, n_targets_sig_down = NA_integer_,
                  mean_log2FC = NA_real_, median_log2FC = NA_real_, max_log2FC = NA_real_, score = NA_real_))
  }
  tibble(
    n_targets_detected = n,
    n_targets_sig_up = sum(sub$padj < SIG_PADJ & sub$log2FC > SIG_LOG2FC, na.rm = TRUE),
    n_targets_sig_down = sum(sub$padj < SIG_PADJ & sub$log2FC < -SIG_LOG2FC, na.rm = TRUE),
    mean_log2FC = mean(sub$log2FC, na.rm = TRUE),
    median_log2FC = median(sub$log2FC, na.rm = TRUE),
    max_log2FC = max(sub$log2FC, na.rm = TRUE),
    score = mean(sub$log2FC, na.rm = TRUE) * log2(n + 1)
  )
}
