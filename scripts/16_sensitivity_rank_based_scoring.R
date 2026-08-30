###############################################################################
# Sensitivity analysis: rank-based cross-dataset integration
#
# Purpose:
#   This analysis provides (a) the exact scoring rule used to build the
#   candidate shortlist (documented in Methods 2.3 / candidate_biomarker
#   _table.R) and (b) a demonstration that the resulting top candidates are
#   not an artefact of that particular scoring rule, via an independent,
#   study-level signed/rank-based integration that does not treat every
#   contrast as an equivalent evidence unit.
#
#   This script builds an independent, rank-based composite score:
#     - For EACH contrast, and among ALL genome-wide detected genes in that
#       contrast (not just the curated ~57-gene signature panel), compute the
#       percentile rank of log2FC (0 = most downregulated, 1 = most
#       upregulated). This normalizes each contrast to its own effect-size
#       distribution before combining across contrasts, which is exactly the
#       "do not treat every contrast as an equivalent evidence unit" goal
#       above (a HUVEC microarray log2FC and an AGM RNA-seq
#       log2FC are not on the same natural scale; percentile rank puts them
#       on a common, comparable footing).
#     - Multiple probes mapping to the same gene symbol within a contrast are
#       collapsed exactly as in the primary pipeline: smallest adjusted P
#       value retained, ties broken by largest |log2FC| (see Methods 2.3).
#     - The GSE33133 duplicate HUVEC contrast (GSE33133_HUVEC_NiV_vs_Mock) is
#       excluded from scoring for the same reason applied throughout this project (Methods 2.3).
#     - A gene's rank_based_score is the mean percentile rank across the
#       scored contrasts in which it was detected (genes not detected in a
#       contrast are excluded from that contrast's average, consistent with
#       the "not detected" missing-data handling already described in
#       Methods 2.3).
#
#   This rank-based score is a SECONDARY robustness check only. It does not
#   replace the primary priority-score-based shortlist (Table S5) used
#   throughout the Results/Figures/Tables; it is reported here solely to
#   demonstrate that the shortlisted candidates are not scoring-rule-specific.
#
# Inputs:
#   GSE32902/results/tables/GSE32902_limma_all_probes_annotated.csv
#   GSE33133/results/tables/GSE33133_limma_all_{NiV_vs_Mock,NiVdC_vs_Mock,NiVdC_vs_NiV}_annotated.csv
#   GSE310471/results/tables/GSE310471_DESeq2_all_{Lung,Tonsil}_{3,4,5}DPI_vs_baseline.csv
#   integrated_results/tables/integrated_signature_long.csv   (curated gene/category panel + scoring_included flag)
#   advanced_analyses/tables/candidate_biomarker_table_full.csv
#   advanced_analyses/tables/candidate_biomarker_shortlist.csv
#
# Outputs:
#   advanced_analyses/tables/sensitivity_rank_based_scoring_full.csv
#   advanced_analyses/tables/sensitivity_shortlist_overlap_summary.csv
#   advanced_analyses/figures/sensitivity_priority_vs_rank_based_scatter.png
#   Console summary (Spearman correlation, per-module shortlist overlap)
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

###############################################################################
# 1. Genome-wide per-contrast tables (same six comparisons as the primary
#    integration script), read WITHOUT restricting to the curated gene panel.
###############################################################################

gse32902_table <- file.path(project_dir, "GSE32902/results/tables/GSE32902_limma_all_probes_annotated.csv")

gse33133_tables <- c(
  GSE33133_HUVEC_NiV_vs_Mock =
    file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiV_vs_Mock_annotated.csv"),
  GSE33133_HUVEC_NiVdC_vs_Mock =
    file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_Mock_annotated.csv"),
  GSE33133_HUVEC_NiVdC_vs_NiV =
    file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_NiV_annotated.csv")
)

gse310471_tables <- c(
  GSE310471_Lung_3DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_3DPI_vs_baseline.csv"),
  GSE310471_Lung_4DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_4DPI_vs_baseline.csv"),
  GSE310471_Lung_5DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_5DPI_vs_baseline.csv"),
  GSE310471_Tonsil_3DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_3DPI_vs_baseline.csv"),
  GSE310471_Tonsil_4DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_4DPI_vs_baseline.csv"),
  GSE310471_Tonsil_5DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_5DPI_vs_baseline.csv")
)

# Same duplicate-contrast exclusion applied throughout this project (Methods 2.3).
scoring_excluded_contrasts <- c("GSE33133_HUVEC_NiV_vs_Mock")

collapse_probes <- function(df, gene_col, log2fc_col, padj_col) {
  df %>%
    mutate(
      gene = str_trim(as.character(.data[[gene_col]])),
      log2FC = as.numeric(.data[[log2fc_col]]),
      padj = as.numeric(.data[[padj_col]])
    ) %>%
    mutate(gene = na_if(gene, "")) %>%
    filter(!is.na(gene), !is.na(log2FC), !is.na(padj)) %>%
    arrange(gene, padj, dplyr::desc(abs(log2FC))) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene, log2FC, padj)
}

genome_wide_ranks <- list()

# GSE32902
d32902 <- readr::read_csv(gse32902_table, show_col_types = FALSE) %>%
  collapse_probes("SYMBOL", "logFC", "adj.P.Val") %>%
  mutate(contrast = "GSE32902_HUVEC_NiV_vs_Mock")
genome_wide_ranks[["GSE32902_HUVEC_NiV_vs_Mock"]] <- d32902

# GSE33133 (3 contrasts)
for (cn in names(gse33133_tables)) {
  d <- readr::read_csv(gse33133_tables[[cn]], show_col_types = FALSE) %>%
    collapse_probes("SYMBOL", "logFC", "adj.P.Val") %>%
    mutate(contrast = cn)
  genome_wide_ranks[[cn]] <- d
}

# GSE310471 (6 contrasts)
for (cn in names(gse310471_tables)) {
  d <- readr::read_csv(gse310471_tables[[cn]], show_col_types = FALSE) %>%
    collapse_probes("Gene", "log2FoldChange", "padj") %>%
    mutate(contrast = cn)
  genome_wide_ranks[[cn]] <- d
}

genome_wide_long <- bind_rows(genome_wide_ranks) %>%
  filter(!contrast %in% scoring_excluded_contrasts) %>%
  group_by(contrast) %>%
  mutate(
    n_detected_genome_wide = n(),
    pct_rank_signed = percent_rank(log2FC)   # 0 = most downregulated, 1 = most upregulated, genome-wide within this contrast
  ) %>%
  ungroup()

readr::write_csv(
  genome_wide_long %>% select(gene, contrast, log2FC, padj, pct_rank_signed, n_detected_genome_wide),
  file.path(table_dir, "sensitivity_genome_wide_percentile_ranks_long.csv")
)

###############################################################################
# 2. Restrict to the curated signature panel already reported in Table S4,
#    and aggregate into one rank-based composite score per gene.
###############################################################################

signature_long <- readr::read_csv(
  file.path(integrated_dir, "tables", "integrated_signature_long.csv"),
  show_col_types = FALSE
)

curated_genes <- signature_long %>% distinct(gene, category)

rank_based_scores <- genome_wide_long %>%
  inner_join(curated_genes, by = "gene") %>%
  group_by(gene, category) %>%
  summarise(
    n_contrasts_detected = n(),
    rank_based_score = mean(pct_rank_signed),
    .groups = "drop"
  )

###############################################################################
# 3. Compare against the primary priority-score shortlist.
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
  left_join(rank_based_scores %>% select(-category), by = "gene") %>%
  mutate(
    priority_rank = rank(dplyr::desc(priority_score), ties.method = "min"),
    rank_based_rank = rank(dplyr::desc(rank_based_score), ties.method = "min"),
    in_priority_shortlist = gene %in% candidate_shortlist$gene
  ) %>%
  arrange(priority_rank)

readr::write_csv(comparison, file.path(table_dir, "sensitivity_rank_based_scoring_full.csv"))

overall_spearman <- cor(
  comparison$priority_score, comparison$rank_based_score,
  method = "spearman", use = "complete.obs"
)

# Per-module: does the rank-based score recover the same shortlist genes,
# using the SAME shortlist size per module as the primary pipeline (n = 15,
# restricted to the same eligible manuscript_module / evidence_strength pool
# candidate_biomarker_table.R already applies)?
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
  left_join(rank_based_scores, by = "gene")

overlap_summary <- eligible_pool %>%
  group_by(manuscript_module) %>%
  group_modify(function(df, key) {
    shortlist_genes <- candidate_shortlist$gene[candidate_shortlist$manuscript_module == key$manuscript_module]
    top_rank_based <- df %>%
      filter(!is.na(rank_based_score)) %>%
      slice_max(rank_based_score, n = length(shortlist_genes), with_ties = FALSE) %>%
      pull(gene)
    n_overlap <- length(intersect(shortlist_genes, top_rank_based))
    n_union <- length(union(shortlist_genes, top_rank_based))
    tibble(
      n_priority_shortlist = length(shortlist_genes),
      n_rank_based_top = length(top_rank_based),
      n_overlap = n_overlap,
      jaccard = ifelse(n_union > 0, n_overlap / n_union, NA_real_),
      priority_only = paste(setdiff(shortlist_genes, top_rank_based), collapse = "; "),
      rank_based_only = paste(setdiff(top_rank_based, shortlist_genes), collapse = "; ")
    )
  }) %>%
  ungroup()

readr::write_csv(overlap_summary, file.path(table_dir, "sensitivity_shortlist_overlap_summary.csv"))

###############################################################################
# 4. Scatter plot: priority score vs rank-based score
###############################################################################

plot_df <- comparison %>%
  filter(!is.na(rank_based_score)) %>%
  mutate(
    label_gene = ifelse(in_priority_shortlist, gene, NA_character_)
  )

p <- ggplot(plot_df, aes(x = rank_based_score, y = priority_score, color = manuscript_module)) +
  geom_point(size = 2, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = label_gene), size = 2.6, max.overlaps = 40, show.legend = FALSE) +
  labs(
    title = "Sensitivity check: priority score vs. rank-based composite score",
    subtitle = paste0("Spearman rho = ", round(overall_spearman, 2), " (labeled points = Table S5 shortlist genes)"),
    x = "Rank-based composite score (mean genome-wide percentile rank of log2FC across detected, scored contrasts)",
    y = "Priority score (primary scoring rule, Methods 2.3)",
    color = "Manuscript module"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))

ggsave(
  file.path(figure_dir, "sensitivity_priority_vs_rank_based_scatter.png"),
  p, width = 9.5, height = 7.5, dpi = 300
)

###############################################################################
# 5. Console summary
###############################################################################

cat("\n=== Sensitivity analysis: rank-based vs. priority score ===\n")
cat("Overall Spearman rho (priority_score vs rank_based_score):", round(overall_spearman, 3), "\n\n")
cat("Per-module shortlist overlap (priority-score shortlist vs. top-N by rank-based score):\n")
print(as.data.frame(overlap_summary))
cat("\nFull gene-level comparison written to: sensitivity_rank_based_scoring_full.csv\n")
cat("Overlap summary written to: sensitivity_shortlist_overlap_summary.csv\n")
cat("Scatter plot written to: sensitivity_priority_vs_rank_based_scatter.png\n")
