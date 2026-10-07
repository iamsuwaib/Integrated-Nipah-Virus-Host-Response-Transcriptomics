###############################################################################
# Cell-composition marker-gene assessment
#
# Purpose:
#   The HUVEC-vs-in-vivo comparison changes cell type, species, platform, and
#   infection timing simultaneously, so bulk tissue complement/coagulation
#   signatures could reflect altered leukocyte, endothelial, or stromal
#   abundance rather than (or in addition to) genuine within-cell
#   transcriptional regulation. This script performs a systematic
#   marker-based assessment of endothelial, myeloid, lymphoid, stromal, pan-immune and epithelial
#   populations as a lighter-weight alternative to formal deconvolution.
#
#   This script computes, per GSE310471 lung/tonsil sample, a composition
#   score for each of six canonical marker-gene sets, then asks two
#   questions per tissue:
#     (1) Does the marker score track infection time (DPI)? (Spearman
#         correlation vs numeric DPI: baseline = 0.)
#     (2) Does the marker score track the tissue's disease/antiviral WGCNA
#         module eigengene reported in the manuscript (lung tan/ME12
#         antiviral module; tonsil blue/ME2 complement/coagulation module)?
#         A strong correlation here would indicate the module's signal is at
#         least partly attributable to compositional shift in that marker
#         population, rather than purely within-cell regulation; a weak
#         correlation does not rule out composition as a contributor (bulk
#         RNA-seq cannot fully separate the two) but is at least consistent
#         with a within-cell regulatory component.
#
#   This is reported as a hypothesis-generating, marker-based screen, not a
#   formal deconvolution, and does not replace or override any WGCNA,
#   integration, or prioritization result reported elsewhere in the study.
#
# Scoring and checks:
#   Genes differ greatly in baseline abundance, so on the raw log scale a highly
#   expressed gene can dominate the mean of a small marker set. This script
#   therefore:
#     (a) standardizes each gene (z-score across that tissue's samples, on
#         the log2 scale) BEFORE averaging within a marker set, so every
#         gene contributes on a comparable scale regardless of its baseline
#         expression level; this is the primary composition score;
#     (b) additionally computes a rank-based single-sample score for every
#         marker set, as a secondary sensitivity cross-check against the
#         z-standardized score;
#     (c) applies Benjamini-Hochberg FDR correction across the marker-set
#         correlation tests (separately for the DPI family and the disease-
#         module family of tests), reported alongside the nominal P values;
#     (d) examines a panel of six canonical lung-specific markers (SFTPC,
#         SFTPB, SFTPA1, NAPSA, SCGB1A1, AGER) plus the broad epithelial
#         marker EPCAM, to test whether an infection-associated, lung-specific
#         signal in tonsil is isolated to SFTPC or reflects broader lung-derived
#         transcript content.
#
# Marker sets (well-established canonical markers present in the GSE310471
# count matrix):
#   endothelial: PECAM1, CDH5, VWF, CLDN5
#   myeloid (monocyte/macrophage/DC): CD68, ITGAM, CD14, LYZ
#   lymphoid (T/B/NK): CD3D, CD3E, MS4A1, NKG7, CD8A
#   pan_immune (leukocyte): PTPRC
#   stromal (fibroblast): COL1A1, COL1A2, DCN, LUM
#   epithelial: EPCAM, SFTPC (SFTPC is lung-alveolar-specific; expected to be
#     uninformative/low in tonsil, which is included for completeness rather
#     than because it is anatomically expected to be meaningful there - see
#     the broader lung-specific-marker investigation below)
#
# Broader lung-specific marker panel (follow-up investigation only, not part
# of the "epithelial" marker set used in the main composition-score table):
#   SFTPC  (alveolar type-II pneumocyte, surfactant protein C)
#   SFTPB  (alveolar type-II pneumocyte, surfactant protein B)
#   SFTPA1 (alveolar type-II pneumocyte, surfactant protein A1)
#   NAPSA  (alveolar type-II pneumocyte, napsin A)
#   SCGB1A1 (airway club/secretory cell, secretoglobin 1A1)
#   AGER   (alveolar type-I pneumocyte, advanced glycosylation end-product receptor)
#
# Inputs:
#   GSE310471/results/tables/GSE310471_Lung_normalized_counts.csv
#   GSE310471/results/tables/GSE310471_Tonsil_normalized_counts.csv
#   GSE310471/results/tables/GSE310471_sample_metadata.csv
#   advanced_analyses/tables/WGCNA_Lung_object.rds   (for ME12, tan antiviral module)
#   advanced_analyses/tables/WGCNA_Tonsil_object.rds (for ME2, blue complement/coagulation module)
#
# Outputs:
#   advanced_analyses/tables/cell_composition_marker_scores_per_sample.csv
#   advanced_analyses/tables/cell_composition_marker_correlation_summary.csv
#   advanced_analyses/tables/cell_composition_marker_correlation_summary_rank_based_sensitivity.csv
#   advanced_analyses/figures/cell_composition_marker_scores_by_dpi.png
#   advanced_analyses/tables/cell_composition_lung_specific_marker_breakdown.csv
#   advanced_analyses/tables/cell_composition_lung_specific_marker_correlation.csv
#   advanced_analyses/figures/cell_composition_lung_specific_marker_breakdown.png
#   Console summary
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "WGCNA", "ggplot2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "WGCNA") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

library(tidyverse)
library(WGCNA)
library(ggplot2)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
gse_dir <- file.path(project_dir, "GSE310471", "results", "tables")
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

marker_sets <- list(
  endothelial = c("PECAM1", "CDH5", "VWF", "CLDN5"),
  myeloid     = c("CD68", "ITGAM", "CD14", "LYZ"),
  lymphoid    = c("CD3D", "CD3E", "MS4A1", "NKG7", "CD8A"),
  pan_immune  = c("PTPRC"),
  stromal     = c("COL1A1", "COL1A2", "DCN", "LUM"),
  epithelial  = c("EPCAM", "SFTPC")
)

# Broader lung-specific marker panel for the SFTPC follow-up investigation
# Is the tonsil epithelial signal isolated to SFTPC, or does it reflect
# broader lung-to-tonsil cross-contamination?
lung_specific_markers_extended <- c("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "SCGB1A1", "AGER")

metadata <- readr::read_csv(file.path(gse_dir, "GSE310471_sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(
    dpi_stripped = ifelse(dpi == "baseline", "0", str_remove(dpi, "DPI")),
    dpi_numeric = as.numeric(dpi_stripped)
  ) %>%
  select(-dpi_stripped)

# Gene-wise z-standardization across samples:
# each gene's log2(count + 1) values are centered/scaled across the samples
# of ONE tissue before any set-level averaging, so no single highly-expressed
# gene can dominate a marker set's score purely due to its baseline level.
zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  dimnames(z) <- dimnames(mat)
  z
}

# Rank-based single-sample score: for each gene,
# rescale the within-tissue sample ranks of log2(count + 1) to [0, 1].
# Reported as a secondary sensitivity cross-check against the z-standardized
# score.
rank01_rows <- function(mat) {
  r <- t(apply(mat, 1, function(x) {
    if (length(unique(x)) == 1) return(rep(0.5, length(x)))
    (rank(x, ties.method = "average") - 1) / (length(x) - 1)
  }))
  dimnames(r) <- dimnames(mat)
  r
}

score_tissue <- function(tissue_name, counts_path) {
  counts <- readr::read_csv(counts_path, show_col_types = FALSE)

  count_mat <- counts %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  log_counts <- log2(count_mat + 1)
  z_counts <- zscore_rows(log_counts)
  rank_counts <- rank01_rows(log_counts)

  scores <- map_dfr(names(marker_sets), function(set_name) {
    genes <- marker_sets[[set_name]]
    genes_present <- intersect(genes, rownames(log_counts))
    if (length(genes_present) == 0) return(NULL)

    tibble(
      sample_id = colnames(log_counts),
      marker_set = set_name,
      n_genes_used = length(genes_present),
      score = as.numeric(colMeans(z_counts[genes_present, , drop = FALSE])),
      score_rank_based = as.numeric(colMeans(rank_counts[genes_present, , drop = FALSE])),
      score_raw_mean_log2 = as.numeric(colMeans(log_counts[genes_present, , drop = FALSE]))
    )
  })

  scores %>% mutate(tissue = tissue_name)
}

lung_scores <- score_tissue("Lung", file.path(gse_dir, "GSE310471_Lung_normalized_counts.csv"))
tonsil_scores <- score_tissue("Tonsil", file.path(gse_dir, "GSE310471_Tonsil_normalized_counts.csv"))

all_scores <- bind_rows(lung_scores, tonsil_scores) %>%
  left_join(metadata %>% select(sample_id, dpi, dpi_numeric, group), by = "sample_id")

###############################################################################
# Disease/antiviral module eigengenes: lung tan/ME12 (antiviral), tonsil
# blue/ME2 (complement/coagulation), as identified in the manuscript's
# WGCNA analysis (advanced_analyses/tables/WGCNA_{Lung,Tonsil}_gene_modules.csv).
###############################################################################

load_module_eigengene <- function(rds_path, me_col) {
  obj <- readRDS(rds_path)
  MEs <- obj$MEs
  if (!me_col %in% colnames(MEs)) {
    stop("Module eigengene column '", me_col, "' not found in ", rds_path,
         ". Available columns: ", paste(colnames(MEs), collapse = ", "))
  }
  tibble(
    sample_id = rownames(MEs),
    module_eigengene = MEs[[me_col]]
  )
}

lung_me <- load_module_eigengene(file.path(table_dir, "WGCNA_Lung_object.rds"), "ME12") %>%
  mutate(tissue = "Lung", module_label = "tan/ME12 (lung antiviral module)")
tonsil_me <- load_module_eigengene(file.path(table_dir, "WGCNA_Tonsil_object.rds"), "ME2") %>%
  mutate(tissue = "Tonsil", module_label = "blue/ME2 (tonsil complement/coagulation module)")

module_eigengenes <- bind_rows(lung_me, tonsil_me)

all_scores <- all_scores %>% left_join(module_eigengenes %>% select(sample_id, module_eigengene, module_label), by = "sample_id")

readr::write_csv(all_scores, file.path(table_dir, "cell_composition_marker_scores_per_sample.csv"))

###############################################################################
# Correlation summary (primary: z-standardized score), with BH/FDR correction
# applied across the marker-set x tissue tests.
# A parallel summary using the rank-based score is written separately as a
# sensitivity cross-check.
###############################################################################

summarise_correlations <- function(df, score_col) {
  df %>%
    group_by(tissue, marker_set) %>%
    summarise(
      n_samples = n(),
      module_label = dplyr::first(module_label),
      rho_vs_DPI = as.numeric(stats::cor(as.numeric(.data[[score_col]]), as.numeric(dpi_numeric), method = "spearman")),
      p_vs_DPI = as.numeric(stats::cor.test(as.numeric(.data[[score_col]]), as.numeric(dpi_numeric), method = "spearman", exact = FALSE)$p.value),
      rho_vs_disease_module = as.numeric(stats::cor(as.numeric(.data[[score_col]]), as.numeric(module_eigengene), method = "spearman")),
      p_vs_disease_module = as.numeric(stats::cor.test(as.numeric(.data[[score_col]]), as.numeric(module_eigengene), method = "spearman", exact = FALSE)$p.value),
      .groups = "drop"
    ) %>%
    arrange(tissue, marker_set) %>%
    mutate(
      # Each family of tests (all DPI tests; all disease-module tests) is
      # corrected separately, since they answer two distinct questions.
      p_vs_DPI_fdr = p.adjust(p_vs_DPI, method = "BH"),
      p_vs_disease_module_fdr = p.adjust(p_vs_disease_module, method = "BH")
    )
}

correlation_summary <- summarise_correlations(all_scores, "score")
readr::write_csv(correlation_summary, file.path(table_dir, "cell_composition_marker_correlation_summary.csv"))

correlation_summary_rank_based <- summarise_correlations(all_scores, "score_rank_based")
readr::write_csv(
  correlation_summary_rank_based,
  file.path(table_dir, "cell_composition_marker_correlation_summary_rank_based_sensitivity.csv")
)

score_vs_rank_agreement <- all_scores %>%
  group_by(tissue, marker_set) %>%
  summarise(
    rho_score_vs_rank_score = as.numeric(stats::cor(score, score_rank_based, method = "spearman")),
    .groups = "drop"
  )

###############################################################################
# Plot: marker scores by DPI, faceted by tissue and marker set
# (plotting the primary z-standardized score)
###############################################################################

plot_df <- all_scores %>%
  mutate(dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI")))

p <- ggplot(plot_df, aes(x = dpi, y = score, color = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.8) +
  facet_grid(marker_set ~ tissue, scales = "free_y") +
  labs(
    title = "Cell-composition marker-gene scores across infection time",
    subtitle = "Mean within-tissue gene-wise z-score (log2 scale) across each canonical marker set, by tissue and DPI",
    x = "Days post infection",
    y = "Marker-set score (mean gene-wise z-score)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0))

ggsave(
  file.path(figure_dir, "cell_composition_marker_scores_by_dpi.png"),
  p, width = 8, height = 11, dpi = 300
)

###############################################################################
# Broader lung-specific marker investigation
#
# The "epithelial" marker set (EPCAM, SFTPC) showed an unexpectedly
# strong correlation with both DPI and the tonsil disease module. SFTPC is a
# lung alveolar type-II pneumocyte marker with no expected role in tonsil, so
# this section tests whether that pattern is specific to SFTPC alone or
# reflects a broader lung-to-tonsil mapping/counting artifact, by adding five
# further canonical lung-specific markers (SFTPB, SFTPA1, NAPSA, SCGB1A1,
# AGER) alongside SFTPC and the broad epithelial marker EPCAM.
###############################################################################

genes_of_interest <- unique(c("EPCAM", lung_specific_markers_extended))

lung_specific_breakdown <- bind_rows(
  lapply(c("Lung", "Tonsil"), function(tissue_name) {
    counts_path <- if (tissue_name == "Lung") {
      file.path(gse_dir, "GSE310471_Lung_normalized_counts.csv")
    } else {
      file.path(gse_dir, "GSE310471_Tonsil_normalized_counts.csv")
    }
    counts <- readr::read_csv(counts_path, show_col_types = FALSE)
    log_counts <- counts %>% column_to_rownames("Gene") %>% as.matrix()
    log_counts <- log2(log_counts + 1)
    z_counts <- zscore_rows(log_counts)

    genes_present <- intersect(genes_of_interest, rownames(log_counts))
    map_dfr(genes_present, function(g) {
      tibble(
        sample_id = colnames(log_counts),
        gene = g,
        log2_count = as.numeric(log_counts[g, ]),
        z_score = as.numeric(z_counts[g, ])
      )
    }) %>% mutate(tissue = tissue_name)
  })
) %>%
  left_join(metadata %>% select(sample_id, dpi, dpi_numeric), by = "sample_id") %>%
  left_join(module_eigengenes %>% select(sample_id, module_eigengene, module_label), by = "sample_id") %>%
  mutate(
    marker_category = ifelse(gene == "EPCAM", "broad epithelial", "lung-specific")
  )

readr::write_csv(lung_specific_breakdown, file.path(table_dir, "cell_composition_lung_specific_marker_breakdown.csv"))

lung_specific_correlation <- lung_specific_breakdown %>%
  group_by(tissue, gene, marker_category) %>%
  summarise(
    n_samples = n(),
    module_label = dplyr::first(module_label),
    rho_vs_DPI = as.numeric(stats::cor(as.numeric(z_score), as.numeric(dpi_numeric), method = "spearman")),
    p_vs_DPI = as.numeric(stats::cor.test(as.numeric(z_score), as.numeric(dpi_numeric), method = "spearman", exact = FALSE)$p.value),
    rho_vs_disease_module = as.numeric(stats::cor(as.numeric(z_score), as.numeric(module_eigengene), method = "spearman")),
    p_vs_disease_module = as.numeric(stats::cor.test(as.numeric(z_score), as.numeric(module_eigengene), method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  arrange(tissue, gene) %>%
  mutate(
    p_vs_DPI_fdr = p.adjust(p_vs_DPI, method = "BH"),
    p_vs_disease_module_fdr = p.adjust(p_vs_disease_module, method = "BH")
  )

readr::write_csv(lung_specific_correlation, file.path(table_dir, "cell_composition_lung_specific_marker_correlation.csv"))

p2 <- ggplot(
  lung_specific_breakdown %>% mutate(dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI"))),
  aes(x = dpi, y = z_score, color = tissue)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.8) +
  facet_grid(gene ~ tissue, scales = "free_y") +
  labs(
    title = "Broad epithelial vs. lung-specific markers, by tissue",
    subtitle = "EPCAM (broad epithelial) vs. six canonical lung-specific markers (within-tissue gene-wise z-score)",
    x = "Days post infection",
    y = "Gene-wise z-score (log2 scale)"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0))

ggsave(
  file.path(figure_dir, "cell_composition_lung_specific_marker_breakdown.png"),
  p2, width = 7, height = max(6, 1.1 * length(genes_of_interest)), dpi = 300, limitsize = FALSE
)

n_lung_specific_sig_in_tonsil_fdr <- lung_specific_correlation %>%
  filter(tissue == "Tonsil", marker_category == "lung-specific") %>%
  summarise(n_dpi = sum(p_vs_DPI_fdr < 0.05), n_module = sum(p_vs_disease_module_fdr < 0.05), n_total = n())

cat("\n=== Broader lung-specific marker investigation ===\n")
print(as.data.frame(lung_specific_correlation))
cat(
  "\nOf ", n_lung_specific_sig_in_tonsil_fdr$n_total,
  " lung-specific markers tested in tonsil, ", n_lung_specific_sig_in_tonsil_fdr$n_dpi,
  " remain FDR-significant vs. DPI and ", n_lung_specific_sig_in_tonsil_fdr$n_module,
  " vs. the disease module after BH correction.\n",
  "If this count is 1 (SFTPC only), the pattern is consistent with an isolated,\n",
  "gene-specific mapping/counting artifact; if multiple lung-specific markers are\n",
  "jointly implicated, this would instead indicate broader lung-to-tonsil\n",
  "cross-contamination with implications for interpretation of the tonsil\n",
  "transcriptomes (see Discussion).\n",
  sep = ""
)

###############################################################################
# Console summary
###############################################################################

cat("\n=== Cell-composition marker-gene assessment (primary: z-standardized score) ===\n")
print(as.data.frame(correlation_summary))

cat("\n=== Rank-based score (secondary sensitivity cross-check) ===\n")
print(as.data.frame(correlation_summary_rank_based))

cat("\nAgreement between z-standardized and rank-based scores (Spearman rho per tissue x marker set):\n")
print(as.data.frame(score_vs_rank_agreement))

cat("\nPer-sample scores (z-standardized, rank-based, and legacy raw mean) written to: cell_composition_marker_scores_per_sample.csv\n")
cat("Correlation summary (z-standardized, with BH/FDR) written to: cell_composition_marker_correlation_summary.csv\n")
cat("Rank-based sensitivity correlation summary written to: cell_composition_marker_correlation_summary_rank_based_sensitivity.csv\n")
cat("Plot written to: cell_composition_marker_scores_by_dpi.png\n")
cat("Broader lung-specific marker breakdown written to: cell_composition_lung_specific_marker_breakdown.csv\n")
cat("Broader lung-specific marker correlation (with BH/FDR) written to: cell_composition_lung_specific_marker_correlation.csv\n")
cat("Broader lung-specific marker plot written to: cell_composition_lung_specific_marker_breakdown.png\n")
