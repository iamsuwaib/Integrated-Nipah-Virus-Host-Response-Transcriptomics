###############################################################################
# TF / upstream regulator signature scoring
#
# Goal:
#   Score targeted TF/regulatory programs relevant to the current biological
#   story: IRF7, STAT1/STAT2/IRF9, ISGF3-like IFN signaling, NF-kB, and selected
#   antiviral regulators.
#
# This script is intentionally offline/reproducible. It does not depend on IPA,
# Enrichr, Dorothea, or internet access.
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
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

integrated_long <- readr::read_csv(
  file.path(project_dir, "integrated_results/tables/integrated_signature_long.csv"),
  show_col_types = FALSE
)

# scoring_included flags GSE33133_HUVEC_NiV_vs_Mock as a duplicate of GSE32902_HUVEC_NiV_vs_Mock
# (same GEO sample records). Retained for display, excluded from
# the recurrence summary (ranked_regulators) below.
scoring_excluded_contrasts <- integrated_long %>%
  filter(!scoring_included) %>%
  pull(contrast) %>%
  unique()

tf_sets <- tribble(
  ~regulator, ~target_genes,
  "IRF7 / antiviral IRF axis", "MX1,MX2,OAS1,OAS2,OAS3,OASL,IFIT1,IFIT2,IFIT3,IFIT5,IFIH1,DDX58,RSAD2,ISG15,USP18,HERC5,HERC6,PARP9",
  "ISGF3-like STAT1/STAT2/IRF9 axis", "STAT1,STAT2,IRF9,MX1,MX2,OAS1,OAS2,OAS3,OASL,IFIT1,IFIT2,IFIT3,ISG15,USP18,RSAD2",
  "RIG-I/MDA5 sensing axis", "DDX58,IFIH1,IRF7,IRF9,STAT1,STAT2,IFIT1,IFIT2,IFIT3,MX1,OAS1,CXCL10,CXCL11",
  "NF-kB / inflammatory chemokine axis", "NFKB1,RELA,NFKBIA,TNF,IL6,CXCL10,CXCL11,CXCL9,CCL2,CCL5,ICAM1,VCAM1,SELE",
  "Endothelial activation axis", "ICAM1,VCAM1,SELE,ANGPT2,VWF,THBD,SERPINE1,PLAU,PLAUR,F3",
  "Complement/coagulation progression axis", "C3,C4A,C4B,C1QA,C1QB,C1QC,CFB,CFD,CFH,CFI,SERPING1,FGB,FGG,FGA,PROC,PROS1,THBD,SERPINE1"
) %>%
  mutate(target_genes = str_split(target_genes, ",")) %>%
  tidyr::unnest(target_genes) %>%
  mutate(target_genes = str_trim(target_genes))

tf_scores <- integrated_long %>%
  mutate(gene = as.character(gene)) %>%
  inner_join(tf_sets, by = c("gene" = "target_genes")) %>%
  group_by(regulator, contrast, dataset, model, tissue, timepoint, scoring_included) %>%
  summarise(
    n_targets_detected = n_distinct(gene),
    n_targets_sig_up = sum(direction == "up", na.rm = TRUE),
    n_targets_sig_down = sum(direction == "down", na.rm = TRUE),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    median_log2FC = median(log2FC, na.rm = TRUE),
    max_log2FC = max(log2FC, na.rm = TRUE),
    score = mean(log2FC, na.rm = TRUE) * log2(n_targets_detected + 1),
    .groups = "drop"
  ) %>%
  arrange(regulator, contrast)

readr::write_csv(tf_scores, file.path(table_dir, "tf_upstream_regulator_signature_scores.csv"))

score_matrix <- tf_scores %>%
  select(regulator, contrast, score) %>%
  pivot_wider(names_from = contrast, values_from = score) %>%
  column_to_rownames("regulator") %>%
  as.matrix()

# Mark the GSE33133 NiV-vs-Mock column as a shared/duplicate sample set: displayed
# here for transparency but excluded from the ranked_regulators recurrence summary below.
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

message("TF/upstream regulator scoring complete.")
