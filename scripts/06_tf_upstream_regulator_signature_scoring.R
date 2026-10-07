###############################################################################
# TF / upstream regulator signature scoring (targeted regulatory-signature axes)
#
# Goal:
#   Score targeted TF/regulatory programs relevant to the biological story:
#   IRF7, STAT1/STAT2/IRF9, ISGF3-like IFN signaling, RIG-I/MDA5 sensing,
#   NF-kB / inflammatory chemokines, endothelial activation and
#   complement/coagulation.
#
# This script is intentionally offline/reproducible. It does not depend on IPA,
# Enrichr, Dorothea, or internet access.
#
# Scores are computed from the genome-wide per-contrast differential-expression
# tables through the shared helper regulator_axis_common.R (also used by script
# 20), so the observed axis scores here and in script 20 are identical. Axis
# genes need not be members of the 57-gene curated panel; for example the NF-kB
# targets NFKB1, RELA and NFKBIA are not in the panel but are scored whenever
# they are detected in a contrast.
#
# Outputs (advanced_analyses/tables, advanced_analyses/figures):
#   tf_upstream_regulator_signature_scores.csv
#   tf_upstream_regulator_ranked_summary.csv
#   tf_axis_target_detection.csv   (which axis genes were detected per contrast)
#   tf_upstream_regulator_score_heatmap.png
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "pheatmap", "RColorBrewer")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
source(file.path(project_dir, "Revised_manuscript_R2/04_github_repository_files/scripts/regulator_axis_common.R"))

advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Contrast metadata (dataset / model / tissue / timepoint) from the integrated signature table
contrast_meta <- readr::read_csv(
  file.path(project_dir, "integrated_results/tables/integrated_signature_long.csv"),
  show_col_types = FALSE
) %>%
  distinct(contrast, dataset, model, tissue, timepoint)

# All 10 contrasts are scored for display; GSE33133 NiV-vs-Mock is flagged
# scoring_included = FALSE (duplicate of the GSE32902 sample records) and is
# excluded from the ranked summary below.
universes <- load_universes(scored_only = FALSE)

tf_scores <- purrr::map_dfr(axes, function(ax) {
  ax_genes <- tf_sets$target_genes[tf_sets$regulator == ax]
  purrr::map_dfr(names(universes), function(cn) {
    dplyr::bind_cols(tibble(regulator = ax, contrast = cn), score_geneset(universes[[cn]], ax_genes))
  })
}) %>%
  left_join(contrast_meta, by = "contrast") %>%
  left_join(contrast_registry %>% select(contrast, scoring_included), by = "contrast") %>%
  select(regulator, contrast, dataset, model, tissue, timepoint, scoring_included,
         n_targets_detected, n_targets_sig_up, n_targets_sig_down,
         mean_log2FC, median_log2FC, max_log2FC, score) %>%
  arrange(regulator, contrast)

readr::write_csv(tf_scores, file.path(table_dir, "tf_upstream_regulator_signature_scores.csv"))

# Which axis genes were (not) detected in each contrast
target_detection <- purrr::map_dfr(axes, function(ax) {
  ax_genes <- tf_sets$target_genes[tf_sets$regulator == ax]
  purrr::map_dfr(names(universes), function(cn) {
    det <- intersect(ax_genes, universes[[cn]]$gene)
    tibble(
      regulator = ax, contrast = cn,
      n_axis_genes = length(ax_genes), n_detected = length(det),
      genes_not_detected = paste(setdiff(ax_genes, det), collapse = ";")
    )
  })
})
readr::write_csv(target_detection, file.path(table_dir, "tf_axis_target_detection.csv"))

score_matrix <- tf_scores %>%
  select(regulator, contrast, score) %>%
  pivot_wider(names_from = contrast, values_from = score) %>%
  column_to_rownames("regulator") %>%
  as.matrix()

scoring_excluded_contrasts <- contrast_registry$contrast[!contrast_registry$scoring_included]
score_matrix_labels_col <- ifelse(
  colnames(score_matrix) %in% scoring_excluded_contrasts,
  paste0(colnames(score_matrix), "*"),
  colnames(score_matrix)
)

png(
  file.path(figure_dir, "tf_upstream_regulator_score_heatmap.png"),
  width = 4200, height = 2600, res = 300
)
pheatmap(
  score_matrix,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  labels_col = score_matrix_labels_col,
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 11,
  angle_col = 45,
  cellwidth = 105,
  cellheight = 34
)
dev.off()

ranked_regulators <- tf_scores %>%
  filter(scoring_included) %>%
  group_by(regulator) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    n_contrasts_positive = sum(score > 0, na.rm = TRUE),
    n_contrasts_high = sum(
      score > quantile(tf_scores$score[tf_scores$scoring_included], 0.75, na.rm = TRUE),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score), desc(n_contrasts_positive))

readr::write_csv(ranked_regulators, file.path(table_dir, "tf_upstream_regulator_ranked_summary.csv"))

message("TF/upstream regulator scoring complete (genome-wide implementation).")
cat("\n=== Ranked regulatory-signature axes (9 scored contrasts) ===\n")
print(as.data.frame(ranked_regulators), row.names = FALSE)
cat("\nNF-kB authoritative mean score:",
    round(ranked_regulators$mean_score[grepl("^NF-kB", ranked_regulators$regulator)], 3), "\n")
