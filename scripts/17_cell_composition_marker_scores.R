###############################################################################
# Cell-composition marker-gene assessment
#
# Purpose:
#   The HUVEC-vs-in-vivo comparison changes cell type, species, platform, and
#   infection timing simultaneously, so bulk tissue complement/coagulation
#   signatures could reflect altered leukocyte, endothelial, or stromal
#   abundance rather than (or in addition to) genuine within-cell
#   transcriptional regulation. This script performs a systematic
#   marker-based assessment of endothelial, myeloid, lymphoid, and stromal
#   populations as a lighter-weight alternative to formal deconvolution.
#
#   This script computes, per GSE310471 lung/tonsil sample, a simple
#   composition score for each of five canonical marker-gene sets (mean
#   log2(normalized count + 1) across the set's genes), then asks two
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
# Marker sets (well-established canonical markers present in the GSE310471
# count matrix):
#   endothelial: PECAM1, CDH5, VWF, CLDN5
#   myeloid (monocyte/macrophage/DC): CD68, ITGAM, CD14, LYZ
#   lymphoid (T/B/NK): CD3D, CD3E, MS4A1, NKG7, CD8A
#   pan_immune (leukocyte): PTPRC
#   stromal (fibroblast): COL1A1, COL1A2, DCN, LUM
#   epithelial: EPCAM, SFTPC (SFTPC is lung-alveolar-specific; expected to be
#     uninformative/low in tonsil, which is included for completeness rather
#     than because it is anatomically expected to be meaningful there)
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
#   advanced_analyses/figures/cell_composition_marker_scores_by_dpi.png
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

metadata <- readr::read_csv(file.path(gse_dir, "GSE310471_sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(
    dpi_stripped = ifelse(dpi == "baseline", "0", str_remove(dpi, "DPI")),
    dpi_numeric = as.numeric(dpi_stripped)
  ) %>%
  select(-dpi_stripped)

score_tissue <- function(tissue_name, counts_path) {
  counts <- readr::read_csv(counts_path, show_col_types = FALSE)
  sample_cols <- setdiff(colnames(counts), "Gene")

  count_mat <- counts %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  log_counts <- log2(count_mat + 1)

  scores <- map_dfr(names(marker_sets), function(set_name) {
    genes <- marker_sets[[set_name]]
    genes_present <- intersect(genes, rownames(log_counts))
    if (length(genes_present) == 0) return(NULL)
    set_scores <- colMeans(log_counts[genes_present, , drop = FALSE])
    tibble(
      sample_id = names(set_scores),
      marker_set = set_name,
      n_genes_used = length(genes_present),
      score = as.numeric(set_scores)
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
# Correlation summary
###############################################################################

correlation_summary <- all_scores %>%
  group_by(tissue, marker_set) %>%
  summarise(
    n_samples = n(),
    module_label = dplyr::first(module_label),
    rho_vs_DPI = as.numeric(stats::cor(as.numeric(score), as.numeric(dpi_numeric), method = "spearman")),
    p_vs_DPI = as.numeric(stats::cor.test(as.numeric(score), as.numeric(dpi_numeric), method = "spearman", exact = FALSE)$p.value),
    rho_vs_disease_module = as.numeric(stats::cor(as.numeric(score), as.numeric(module_eigengene), method = "spearman")),
    p_vs_disease_module = as.numeric(stats::cor.test(as.numeric(score), as.numeric(module_eigengene), method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  arrange(tissue, marker_set)

readr::write_csv(correlation_summary, file.path(table_dir, "cell_composition_marker_correlation_summary.csv"))

###############################################################################
# Plot: marker scores by DPI, faceted by tissue and marker set
###############################################################################

plot_df <- all_scores %>%
  mutate(dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI")))

p <- ggplot(plot_df, aes(x = dpi, y = score, color = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.8) +
  facet_grid(marker_set ~ tissue, scales = "free_y") +
  labs(
    title = "Cell-composition marker-gene scores across infection time",
    subtitle = "Mean log2(normalized count + 1) across each canonical marker set, by tissue and DPI",
    x = "Days post infection",
    y = "Marker-set score"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0))

ggsave(
  file.path(figure_dir, "cell_composition_marker_scores_by_dpi.png"),
  p, width = 8, height = 11, dpi = 300
)

###############################################################################
# Per-gene decomposition of the epithelial marker set (sanity check)
#
# The epithelial set (EPCAM, SFTPC) showed an unexpectedly strong correlation
# with both DPI and the tonsil disease module (see correlation_summary above).
# SFTPC is a lung alveolar type-II pneumocyte marker with no expected role in
# tonsil, so before reporting this result we decompose the two-gene mean score
# into its individual genes to check which one is responsible.
###############################################################################

epithelial_per_gene <- bind_rows(
  lapply(c("Lung", "Tonsil"), function(tissue_name) {
    counts_path <- if (tissue_name == "Lung") {
      file.path(gse_dir, "GSE310471_Lung_normalized_counts.csv")
    } else {
      file.path(gse_dir, "GSE310471_Tonsil_normalized_counts.csv")
    }
    counts <- readr::read_csv(counts_path, show_col_types = FALSE)
    log_counts <- counts %>% column_to_rownames("Gene") %>% as.matrix()
    log_counts <- log2(log_counts + 1)
    genes_present <- intersect(c("EPCAM", "SFTPC"), rownames(log_counts))
    map_dfr(genes_present, function(g) {
      tibble(sample_id = colnames(log_counts), gene = g, log2_count = as.numeric(log_counts[g, ]))
    }) %>% mutate(tissue = tissue_name)
  })
) %>%
  left_join(metadata %>% select(sample_id, dpi, dpi_numeric), by = "sample_id") %>%
  left_join(module_eigengenes %>% select(sample_id, module_eigengene, module_label), by = "sample_id")

readr::write_csv(epithelial_per_gene, file.path(table_dir, "cell_composition_epithelial_per_gene_breakdown.csv"))

epithelial_per_gene_correlation <- epithelial_per_gene %>%
  group_by(tissue, gene) %>%
  summarise(
    n_samples = n(),
    module_label = dplyr::first(module_label),
    rho_vs_DPI = as.numeric(stats::cor(as.numeric(log2_count), as.numeric(dpi_numeric), method = "spearman")),
    p_vs_DPI = as.numeric(stats::cor.test(as.numeric(log2_count), as.numeric(dpi_numeric), method = "spearman", exact = FALSE)$p.value),
    rho_vs_disease_module = as.numeric(stats::cor(as.numeric(log2_count), as.numeric(module_eigengene), method = "spearman")),
    p_vs_disease_module = as.numeric(stats::cor.test(as.numeric(log2_count), as.numeric(module_eigengene), method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  arrange(tissue, gene)

readr::write_csv(epithelial_per_gene_correlation, file.path(table_dir, "cell_composition_epithelial_per_gene_correlation.csv"))

p2 <- ggplot(
  epithelial_per_gene %>% mutate(dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI"))),
  aes(x = dpi, y = log2_count, color = tissue)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.8) +
  facet_grid(gene ~ tissue, scales = "free_y") +
  labs(
    title = "Epithelial marker set decomposed by gene",
    subtitle = "EPCAM (broad epithelial) vs. SFTPC (lung-specific)",
    x = "Days post infection",
    y = "log2(normalized count + 1)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0))

ggsave(
  file.path(figure_dir, "cell_composition_epithelial_per_gene_breakdown.png"),
  p2, width = 7, height = 6, dpi = 300
)

cat("\n=== Epithelial marker set per-gene decomposition (sanity check) ===\n")
print(as.data.frame(epithelial_per_gene_correlation))
cat("\nInterpretation note: if SFTPC (a lung-specific marker) alone accounts for\n")
cat("the tonsil epithelial-set correlation while EPCAM does not, this indicates a\n")
cat("likely technical confound (e.g., cross-tissue material) rather than genuine\n")
cat("tonsil epithelial biology; the epithelial marker set is reported with this\n")
cat("caveat and excluded from biological interpretation (see Discussion).\n")

###############################################################################
# Console summary
###############################################################################

cat("\n=== Cell-composition marker-gene assessment ===\n")
print(as.data.frame(correlation_summary))
cat("\nPer-sample scores written to: cell_composition_marker_scores_per_sample.csv\n")
cat("Correlation summary written to: cell_composition_marker_correlation_summary.csv\n")
cat("Plot written to: cell_composition_marker_scores_by_dpi.png\n")
