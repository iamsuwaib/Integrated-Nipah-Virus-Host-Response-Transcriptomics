###############################################################################
# Permutation-based null model for the targeted regulatory-signature axis
# scores (Response to Reviewers R2-7)
#
# The existing axis scores (02_tf_upstream_regulator_signature_scoring.R,
# score = mean(log2FC) * log2(n_targets_detected + 1)) were reported as raw
# numbers with no significance testing. This script adds a genome-wide
# permutation null: for each contrast, draws many random gene sets of the
# same size as each axis's DETECTED gene count from that contrast's full,
# genome-wide tested-gene universe (NOT from the small curated panel used
# elsewhere for signature integration), recomputes the same score formula,
# and derives an empirical one-sided P value for each axis's observed
# mean_score (averaged across the 9 contrasts used for ranking, matching the
# existing "mean_score" column in the ranked_regulators / Supplementary
# Table S6 summary).
#
# Per-contrast genome-wide gene universes:
#   GSE32902_HUVEC_NiV_vs_Mock        -> GSE32902_limma_all_probes_annotated.csv
#   GSE33133_HUVEC_NiVdC_vs_Mock      -> GSE33133_limma_all_NiVdC_vs_Mock_annotated.csv
#   GSE33133_HUVEC_NiVdC_vs_NiV       -> GSE33133_limma_all_NiVdC_vs_NiV_annotated.csv
#   GSE310471_<Tissue>_<DPI>_vs_baseline -> GSE310471_DESeq2_all_<Tissue>_<DPI>_vs_baseline.csv
#
# GSE33133_HUVEC_NiV_vs_Mock is excluded here exactly as it is excluded from
# ranked_regulators in the main scoring script (duplicate GSE32902 sample
# records; see Response to Reviewers R2-1).
#
# For the microarray (limma) universes, multiple probes can map to the same
# gene SYMBOL. We collapse to one row per SYMBOL using the same rule already
# used elsewhere in this project for the curated panel (integrated_results/
# scripts/integrated_cross_dataset_signature_heatmap.R, read_limma_signature):
# arrange by padj then by descending |log2FC| and keep the first row per
# SYMBOL - applied here genome-wide rather than restricted to the curated
# panel, so every gene gets a fair, non-circular null distribution.
#
# Outputs (advanced_analyses/tables/):
#   regulator_axis_permutation_null_scores.csv - per-axis, per-permutation mean_score (long)
#   regulator_axis_permutation_summary.csv     - per-axis observed mean_score, null mean/SD, empirical P
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
suppressPackageStartupMessages(library(tidyverse))

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")

set.seed(20260830)
N_PERM <- 10000

## --- axis gene sets (identical to 02_tf_upstream_regulator_signature_scoring.R) ---
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

axes <- unique(tf_sets$regulator)

## --- per-contrast genome-wide universes (gene -> log2FC, one row per gene) ---
read_limma_universe <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      SYMBOL = str_trim(SYMBOL),
      SYMBOL = na_if(SYMBOL, "")
    ) %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(SYMBOL, adj.P.Val, desc(abs(logFC))) %>%
    group_by(SYMBOL) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene = SYMBOL, log2FC = as.numeric(logFC))
}

read_deseq2_universe <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(log2FoldChange)) %>%
    transmute(gene = Gene, log2FC = as.numeric(log2FoldChange)) %>%
    distinct(gene, .keep_all = TRUE)
}

contrast_files <- list(
  GSE32902_HUVEC_NiV_vs_Mock = list(
    path = file.path(project_dir, "GSE32902/results/tables/GSE32902_limma_all_probes_annotated.csv"),
    reader = read_limma_universe
  ),
  GSE33133_HUVEC_NiVdC_vs_Mock = list(
    path = file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_Mock_annotated.csv"),
    reader = read_limma_universe
  ),
  GSE33133_HUVEC_NiVdC_vs_NiV = list(
    path = file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_NiV_annotated.csv"),
    reader = read_limma_universe
  ),
  GSE310471_Lung_3DPI_vs_baseline = list(
    path = file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_3DPI_vs_baseline.csv"),
    reader = read_deseq2_universe
  ),
  GSE310471_Lung_4DPI_vs_baseline = list(
    path = file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_4DPI_vs_baseline.csv"),
    reader = read_deseq2_universe
  ),
  GSE310471_Lung_5DPI_vs_baseline = list(
    path = file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_5DPI_vs_baseline.csv"),
    reader = read_deseq2_universe
  ),
  GSE310471_Tonsil_3DPI_vs_baseline = list(
    path = file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_3DPI_vs_baseline.csv"),
    reader = read_deseq2_universe
  ),
  GSE310471_Tonsil_4DPI_vs_baseline = list(
    path = file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_4DPI_vs_baseline.csv"),
    reader = read_deseq2_universe
  ),
  GSE310471_Tonsil_5DPI_vs_baseline = list(
    path = file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_5DPI_vs_baseline.csv"),
    reader = read_deseq2_universe
  )
)

message("Loading per-contrast genome-wide universes...")
universes <- map(contrast_files, function(cf) cf$reader(cf$path))
for (nm in names(universes)) {
  message("  ", nm, ": ", nrow(universes[[nm]]), " genes")
}

score_for_geneset <- function(universe_df, genes) {
  sub <- universe_df %>% filter(gene %in% genes)
  n_detected <- nrow(sub)
  if (n_detected == 0) return(list(score = NA_real_, n_detected = 0))
  list(score = mean(sub$log2FC, na.rm = TRUE) * log2(n_detected + 1), n_detected = n_detected)
}

## --- observed axis scores per contrast (should reproduce 02_...R's tf_scores) ---
observed <- map_dfr(axes, function(ax) {
  ax_genes <- tf_sets$target_genes[tf_sets$regulator == ax]
  map_dfr(names(universes), function(cn) {
    res <- score_for_geneset(universes[[cn]], ax_genes)
    tibble(regulator = ax, contrast = cn, n_detected = res$n_detected, observed_score = res$score)
  })
})
readr::write_csv(observed, file.path(table_dir, "regulator_axis_permutation_observed_scores.csv"))

observed_mean <- observed %>%
  group_by(regulator) %>%
  summarise(observed_mean_score = mean(observed_score, na.rm = TRUE), .groups = "drop")

## --- permutation null: for each replicate, draw a same-size random gene set
##     per contrast (matching n_detected for that axis in that contrast),
##     score it, and average across contrasts ---
message("Running ", N_PERM, " permutations per axis...")

null_rows <- vector("list", length(axes) * N_PERM)
k <- 1
for (ax in axes) {
  ax_obs <- observed %>% filter(regulator == ax)
  message("  axis: ", ax)
  for (p in seq_len(N_PERM)) {
    per_contrast_scores <- map_dbl(names(universes), function(cn) {
      n_det <- ax_obs$n_detected[ax_obs$contrast == cn]
      uni <- universes[[cn]]
      if (length(n_det) == 0 || n_det == 0 || n_det > nrow(uni)) return(NA_real_)
      rand_genes <- uni$log2FC[sample.int(nrow(uni), size = n_det, replace = FALSE)]
      mean(rand_genes, na.rm = TRUE) * log2(n_det + 1)
    })
    null_rows[[k]] <- tibble(
      regulator = ax, perm_id = p,
      null_mean_score = mean(per_contrast_scores, na.rm = TRUE)
    )
    k <- k + 1
  }
}
null_df <- bind_rows(null_rows)
readr::write_csv(null_df, file.path(table_dir, "regulator_axis_permutation_null_scores.csv"))

summary_df <- null_df %>%
  group_by(regulator) %>%
  summarise(
    null_mean = mean(null_mean_score, na.rm = TRUE),
    null_sd = sd(null_mean_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(observed_mean, by = "regulator") %>%
  left_join(
    null_df %>%
      left_join(observed_mean, by = "regulator") %>%
      group_by(regulator) %>%
      summarise(
        empirical_p_one_sided = (sum(null_mean_score >= observed_mean_score, na.rm = TRUE) + 1) / (n() + 1),
        .groups = "drop"
      ),
    by = "regulator"
  ) %>%
  arrange(empirical_p_one_sided)

readr::write_csv(summary_df, file.path(table_dir, "regulator_axis_permutation_summary.csv"))

message("Permutation test complete.")
print(summary_df)
