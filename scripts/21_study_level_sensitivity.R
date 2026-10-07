###############################################################################
# Study-level (equal-weight-by-study) sensitivity analysis + leave-one-
# system-out robustness check
#
# Rationale:
#   The primary priority score (05_candidate_biomarker_table.R) and the
#   existing rank-based sensitivity check (16_sensitivity_rank_based_scoring.R)
#   both treat each of the 9 scored contrasts as an equivalent, independent
#   evidence unit. But the 9 scored contrasts are not 9 independent biological
#   systems: GSE32902_HUVEC_NiV_vs_Mock, GSE33133_HUVEC_NiVdC_vs_Mock and
#   GSE33133_HUVEC_NiVdC_vs_NiV all originate from the same underlying HUVEC
#   experiment (Methods 2.3), and the 3 Lung + 3 Tonsil DPI contrasts all
#   originate from the same GSE310471 animal cohort. A simple mean across
#   contrasts therefore still rewards recurrence within one system rather
#   than independently rewarding breadth across systems.
#
#   This script aggregates evidence WITHIN each of three independent
#   biological systems first - HUVEC (GSE32902 + GSE33133, collapsed into one
#   study), Lung (GSE310471), and Tonsil (GSE310471) - and only then combines
#   the three equally-weighted per-study summaries into a single score. Lung
#   and Tonsil are kept as separate systems (rather than merged into one
#   "GSE310471" unit) because they are distinct tissue compartments and the
#   manuscript's own central claim is that lung- and tonsil-associated disease
#   programs differ; collapsing them would obscure exactly the distinction the
#   paper is making. This accounts for the GSE32902/GSE33133
#   non-independence, since both are merged into one
#   HUVEC study here.
#
#   A formal meta-analysis is not attempted; this is a transparent equal-weight-by-
#   study score plus a leave-one-system-out stability check, reported as a
#   secondary robustness analysis alongside the existing rank-based check
#   (16_sensitivity_rank_based_scoring.R). It does not replace the primary
#   priority-score-based shortlist (Table S5).
#
# Must be run AFTER 16_sensitivity_rank_based_scoring.R, which writes the
# genome-wide per-contrast percentile-rank table this script reads.
#
# Inputs:
#   advanced_analyses/tables/sensitivity_genome_wide_percentile_ranks_long.csv
#     (written by 16_sensitivity_rank_based_scoring.R; already excludes the
#     GSE33133_HUVEC_NiV_vs_Mock duplicate contrast)
#   integrated_results/tables/integrated_signature_long.csv   (curated gene/category panel)
#   advanced_analyses/tables/candidate_biomarker_table_full.csv
#   advanced_analyses/tables/candidate_biomarker_shortlist.csv
#
# Outputs:
#   advanced_analyses/tables/sensitivity_study_level_scoring_full.csv
#   advanced_analyses/tables/sensitivity_study_level_overlap_summary.csv
#   advanced_analyses/figures/sensitivity_study_level_vs_priority_scatter.png
#   advanced_analyses/figures/sensitivity_leave_one_system_out_stability.png
#   Console summary (per-study coverage, Spearman correlations, leave-one-
#   system-out shortlist overlap)
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "ggplot2", "ggrepel")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}

library(tidyverse)
library(ggplot2)
library(ggrepel)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
integrated_dir <- file.path(project_dir, "integrated_results")
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

genome_wide_ranks_path <- file.path(table_dir, "sensitivity_genome_wide_percentile_ranks_long.csv")
if (!file.exists(genome_wide_ranks_path)) {
  stop(
    "Missing ", genome_wide_ranks_path,
    " - run 16_sensitivity_rank_based_scoring.R first (this script reuses its ",
    "genome-wide percentile-rank output rather than recomputing it)."
  )
}

###############################################################################
# 1. Load genome-wide percentile ranks and assign each scored contrast to one
#    of three independent biological systems.
###############################################################################

genome_wide_long <- readr::read_csv(genome_wide_ranks_path, show_col_types = FALSE)

contrast_to_study <- c(
  GSE32902_HUVEC_NiV_vs_Mock        = "HUVEC",
  GSE33133_HUVEC_NiVdC_vs_Mock      = "HUVEC",
  GSE33133_HUVEC_NiVdC_vs_NiV       = "HUVEC",
  GSE310471_Lung_3DPI_vs_baseline   = "Lung",
  GSE310471_Lung_4DPI_vs_baseline   = "Lung",
  GSE310471_Lung_5DPI_vs_baseline   = "Lung",
  GSE310471_Tonsil_3DPI_vs_baseline = "Tonsil",
  GSE310471_Tonsil_4DPI_vs_baseline = "Tonsil",
  GSE310471_Tonsil_5DPI_vs_baseline = "Tonsil"
)

stopifnot(all(unique(genome_wide_long$contrast) %in% names(contrast_to_study)))

genome_wide_long <- genome_wide_long %>%
  mutate(study = unname(contrast_to_study[contrast]))

###############################################################################
# 2. Restrict to the curated signature panel and aggregate WITHIN each study
#    first (mean percentile rank across that study's detected contrasts),
#    then combine the (up to three) per-study means with equal weight.
###############################################################################

signature_long <- readr::read_csv(
  file.path(integrated_dir, "tables", "integrated_signature_long.csv"),
  show_col_types = FALSE
)
curated_genes <- signature_long %>% distinct(gene, category)

per_study_means <- genome_wide_long %>%
  inner_join(curated_genes, by = "gene") %>%
  group_by(gene, category, study) %>%
  summarise(study_mean_pct_rank = mean(pct_rank_signed), .groups = "drop")

study_wide <- per_study_means %>%
  pivot_wider(names_from = study, values_from = study_mean_pct_rank)
for (s in c("HUVEC", "Lung", "Tonsil")) {
  if (!s %in% colnames(study_wide)) study_wide[[s]] <- NA_real_
}

row_mean_ignore_na <- function(...) {
  m <- cbind(...)
  rowMeans(m, na.rm = TRUE)
}

study_level_scores <- study_wide %>%
  mutate(
    n_studies_detected = rowSums(!is.na(cbind(HUVEC, Lung, Tonsil))),
    study_level_score = row_mean_ignore_na(HUVEC, Lung, Tonsil),
    # Leave-one-system-out: equal-weight mean of the REMAINING two studies only.
    loo_drop_HUVEC_score  = row_mean_ignore_na(Lung, Tonsil),
    loo_drop_Lung_score   = row_mean_ignore_na(HUVEC, Tonsil),
    loo_drop_Tonsil_score = row_mean_ignore_na(HUVEC, Lung)
  )

###############################################################################
# 3. Compare against the primary priority-score shortlist (same comparison
#    structure as 16_sensitivity_rank_based_scoring.R).
###############################################################################

candidate_full <- readr::read_csv(
  file.path(table_dir, "candidate_biomarker_table_full.csv"),
  show_col_types = FALSE
)
candidate_shortlist <- readr::read_csv(
  file.path(table_dir, "candidate_biomarker_shortlist.csv"),
  show_col_types = FALSE
)

comparison <- candidate_full %>%
  select(gene, category, manuscript_module, evidence_strength, priority_score) %>%
  left_join(study_level_scores %>% select(-category), by = "gene") %>%
  mutate(
    priority_rank = rank(dplyr::desc(priority_score), ties.method = "min"),
    study_level_rank = rank(dplyr::desc(study_level_score), ties.method = "min"),
    in_priority_shortlist = gene %in% candidate_shortlist$gene
  ) %>%
  arrange(priority_rank)

readr::write_csv(comparison, file.path(table_dir, "sensitivity_study_level_scoring_full.csv"))

score_cols <- c(
  full = "study_level_score",
  drop_HUVEC = "loo_drop_HUVEC_score",
  drop_Lung = "loo_drop_Lung_score",
  drop_Tonsil = "loo_drop_Tonsil_score"
)

eligible_pool <- candidate_full %>%
  filter(
    manuscript_module %in% c(
      "Conserved antiviral/IFN core",
      "In vivo complement/coagulation disease module",
      "HUVEC-enriched endothelial/early-response module"
    ),
    evidence_strength %in% c("high", "moderate_high", "in_vivo_specific", "HUVEC_enriched")
  ) %>%
  select(gene, manuscript_module) %>%
  left_join(study_level_scores, by = "gene")

overlap_summary <- purrr::map_dfr(names(score_cols), function(scenario) {
  score_col <- score_cols[[scenario]]
  eligible_pool %>%
    group_by(manuscript_module) %>%
    group_modify(function(df, key) {
      shortlist_genes <- candidate_shortlist$gene[candidate_shortlist$manuscript_module == key$manuscript_module]
      top_study_level <- df %>%
        filter(!is.na(.data[[score_col]])) %>%
        slice_max(.data[[score_col]], n = length(shortlist_genes), with_ties = FALSE) %>%
        pull(gene)
      n_overlap <- length(intersect(shortlist_genes, top_study_level))
      n_union <- length(union(shortlist_genes, top_study_level))
      tibble(
        scenario = scenario,
        n_priority_shortlist = length(shortlist_genes),
        n_study_level_top = length(top_study_level),
        n_overlap = n_overlap,
        jaccard = ifelse(n_union > 0, n_overlap / n_union, NA_real_),
        priority_only = paste(setdiff(shortlist_genes, top_study_level), collapse = "; "),
        study_level_only = paste(setdiff(top_study_level, shortlist_genes), collapse = "; ")
      )
    }) %>%
    ungroup()
})

readr::write_csv(overlap_summary, file.path(table_dir, "sensitivity_study_level_overlap_summary.csv"))

###############################################################################
# 4. Figure 1: full study-level score vs. priority score (same style as the
#    existing rank-based sensitivity scatter).
###############################################################################

overall_spearman <- cor(
  comparison$priority_score, comparison$study_level_score,
  method = "spearman", use = "complete.obs"
)

plot_df <- comparison %>%
  filter(!is.na(study_level_score)) %>%
  mutate(label_gene = ifelse(in_priority_shortlist, gene, NA_character_))

p1 <- ggplot(plot_df, aes(x = study_level_score, y = priority_score, color = manuscript_module)) +
  geom_point(size = 2, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = label_gene), size = 2.6, max.overlaps = 40, show.legend = FALSE) +
  labs(
    title = "Sensitivity check: priority score vs. equal-weight-by-study composite score",
    subtitle = paste0(
      "Spearman rho = ", round(overall_spearman, 2),
      " (labeled points = Table S5 shortlist genes; score aggregates HUVEC, Lung, Tonsil with equal weight per study)"
    ),
    x = "Study-level composite score (mean of HUVEC / Lung / Tonsil per-study mean percentile rank)",
    y = "Priority score (primary scoring rule, Methods 2.3)",
    color = "Manuscript module"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))

ggsave(
  file.path(figure_dir, "sensitivity_study_level_vs_priority_scatter.png"),
  p1, width = 9.5, height = 7.5, dpi = 300
)

###############################################################################
# 5. Figure 2: leave-one-system-out stability - does each shortlisted gene
#    stay in the top-N (per its module, same N as the primary shortlist) when
#    HUVEC, Lung, or Tonsil is dropped in turn?
###############################################################################

stability_long <- purrr::map_dfr(names(score_cols), function(scenario) {
  score_col <- score_cols[[scenario]]
  eligible_pool %>%
    group_by(manuscript_module) %>%
    group_modify(function(df, key) {
      shortlist_genes <- candidate_shortlist$gene[candidate_shortlist$manuscript_module == key$manuscript_module]
      top_n <- df %>%
        filter(!is.na(.data[[score_col]])) %>%
        slice_max(.data[[score_col]], n = length(shortlist_genes), with_ties = FALSE) %>%
        pull(gene)
      tibble(gene = shortlist_genes, retained = shortlist_genes %in% top_n)
    }) %>%
    ungroup() %>%
    mutate(scenario = scenario)
})

stability_wide <- stability_long %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("full", "drop_HUVEC", "drop_Lung", "drop_Tonsil"),
      labels = c("Full (all 3 systems)", "Drop HUVEC", "Drop Lung", "Drop Tonsil")
    )
  )

gene_order <- candidate_shortlist %>%
  arrange(manuscript_module, desc(priority_score)) %>%
  pull(gene) %>%
  unique()
stability_wide <- stability_wide %>%
  mutate(gene = factor(gene, levels = rev(gene_order)))

p2 <- ggplot(stability_wide, aes(x = scenario, y = gene, fill = retained)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "grey85"),
    labels = c("TRUE" = "In top-N", "FALSE" = "Dropped out"),
    name = NULL
  ) +
  labs(
    title = "Leave-one-system-out stability of the Table S5 shortlist",
    subtitle = "Each column re-ranks genes by the equal-weight-by-study score using only the remaining two systems",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "sensitivity_leave_one_system_out_stability.png"),
  p2, width = 7, height = max(6, 0.11 * length(gene_order)), dpi = 300, limitsize = FALSE
)

###############################################################################
# 6. Console summary
###############################################################################

cat("\n=== Study-level (equal-weight-by-study) sensitivity analysis ===\n")
cat("Per-study gene detection in curated panel:\n")
print(colSums(!is.na(study_level_scores[, c("HUVEC", "Lung", "Tonsil")])))

cat("\nOverall Spearman rho (priority_score vs study_level_score):", round(overall_spearman, 3), "\n")

cat("\nShortlist overlap (priority-score shortlist vs. top-N by study-level score), by scenario:\n")
print(as.data.frame(overlap_summary))

n_total_shortlist <- length(unique(stability_long$gene))
retained_all_three_drops <- stability_long %>%
  filter(scenario != "full") %>%
  group_by(gene) %>%
  summarise(retained_in_all_drops = all(retained), .groups = "drop") %>%
  summarise(n = sum(retained_in_all_drops)) %>%
  pull(n)

cat(
  "\n", retained_all_three_drops, " of ", n_total_shortlist,
  " Table S5 shortlist genes remain in their module's top-N under EVERY ",
  "single-system removal (drop HUVEC / drop Lung / drop Tonsil).\n",
  sep = ""
)

cat("\nFull gene-level comparison written to: sensitivity_study_level_scoring_full.csv\n")
cat("Overlap summary (all 4 scenarios) written to: sensitivity_study_level_overlap_summary.csv\n")
cat("Scatter plot written to: sensitivity_study_level_vs_priority_scatter.png\n")
cat("Leave-one-system-out stability plot written to: sensitivity_leave_one_system_out_stability.png\n")
