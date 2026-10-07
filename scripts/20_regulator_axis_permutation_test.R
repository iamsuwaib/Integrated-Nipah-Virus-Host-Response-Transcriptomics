###############################################################################
# Permutation null model for the targeted regulatory-signature axis scores
#
# Design:
#   1. Gene universes and axis scoring come from the shared helper
#      regulator_axis_common.R (also used by script 06), so the observed axis
#      scores here and in script 06 are identical.
#   2. PRIMARY permutation design = one COMMON random gene set per replicate.
#      In each replicate a single random gene set of the same size as the axis
#      is drawn once from the pooled gene universe (all genes tested in any of
#      the 9 scored contrasts) and then evaluated in all 9 contrasts, using the
#      same formula, with genes counted in a contrast only if they were tested
#      there (exactly as for the observed axis). The replicate statistic is the
#      mean score across the 9 contrasts. This preserves the structure of the
#      question, "does a recurrent biological gene set behave coherently across
#      datasets?", and keeps the correlation between contrasts in the null.
#   3. SENSITIVITY design = an alternative scheme, in which an independent random
#      gene set (matched to the number of axis genes detected in that contrast)
#      is drawn separately for every contrast.
#   4. Multiple testing across the six axes: Benjamini-Hochberg (primary) and
#      Holm (conservative) adjusted empirical P values are reported.
#
# Inputs : genome-wide per-contrast tables (see regulator_axis_common.R)
# Outputs (advanced_analyses/tables):
#   regulator_axis_permutation_observed_scores.csv
#   regulator_axis_permutation_null_scores.csv      (long; both designs)
#   regulator_axis_permutation_summary.csv
###############################################################################

options(stringsAsFactors = FALSE)
packages <- c("tidyverse")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
suppressPackageStartupMessages(library(tidyverse))

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
source(file.path(project_dir, "Revised_manuscript_R2/04_github_repository_files/scripts/regulator_axis_common.R"))
table_dir <- file.path(project_dir, "advanced_analyses", "tables")

set.seed(20261004)
N_PERM <- 10000

message("Loading per-contrast genome-wide universes (9 scored contrasts)...")
universes <- load_universes(scored_only = TRUE)
for (nm in names(universes)) message("  ", nm, ": ", nrow(universes[[nm]]), " genes")

# Pooled universe matrix: rows = genes tested in at least one contrast,
# columns = contrasts, entries = log2FC (NA where the gene was not tested)
all_genes <- sort(unique(unlist(lapply(universes, function(u) u$gene))))
M <- vapply(universes, function(u) u$log2FC[match(all_genes, u$gene)], numeric(length(all_genes)))
rownames(M) <- all_genes
n_contrasts <- ncol(M)
message("Pooled universe: ", nrow(M), " genes x ", n_contrasts, " contrasts")

score_from_matrix <- function(sub) {
  n_det <- colSums(!is.na(sub))
  sc <- ifelse(n_det > 0, colSums(sub, na.rm = TRUE) / pmax(n_det, 1) * log2(n_det + 1), NA_real_)
  list(score = sc, n_det = n_det)
}

###############################################################################
# Observed axis scores
###############################################################################

observed <- purrr::map_dfr(axes, function(ax) {
  ax_genes <- tf_sets$target_genes[tf_sets$regulator == ax]
  purrr::map_dfr(names(universes), function(cn) {
    r <- score_geneset(universes[[cn]], ax_genes)
    tibble(regulator = ax, contrast = cn, n_detected = r$n_targets_detected, observed_score = r$score)
  })
})
readr::write_csv(observed, file.path(table_dir, "regulator_axis_permutation_observed_scores.csv"))

# Internal consistency check: matrix-based scoring must equal score_geneset()
for (ax in axes) {
  ax_genes <- tf_sets$target_genes[tf_sets$regulator == ax]
  sm <- score_from_matrix(M[rownames(M) %in% ax_genes, , drop = FALSE])$score
  so <- observed$observed_score[observed$regulator == ax][match(colnames(M), observed$contrast[observed$regulator == ax])]
  stopifnot(isTRUE(all.equal(unname(sm), unname(so), tolerance = 1e-8)))
}

axis_info <- observed %>%
  group_by(regulator) %>%
  summarise(observed_mean_score = mean(observed_score, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    tf_sets %>% mutate(in_universe = target_genes %in% rownames(M)) %>%
      group_by(regulator) %>%
      summarise(n_axis_genes = n(), n_axis_genes_in_universe = sum(in_universe), .groups = "drop"),
    by = "regulator"
  )

###############################################################################
# Permutation nulls
###############################################################################

null_common_set <- function(K, nperm) {
  out <- numeric(nperm)
  nU <- nrow(M)
  for (p in seq_len(nperm)) {
    sub <- M[sample.int(nU, K), , drop = FALSE]
    out[p] <- mean(score_from_matrix(sub)$score, na.rm = TRUE)
  }
  out
}

# Universe vectors per contrast for the sensitivity (independent per-contrast) design
uni_vec <- lapply(seq_len(n_contrasts), function(j) M[!is.na(M[, j]), j])

null_per_contrast <- function(n_det_vec, nperm) {
  out <- numeric(nperm)
  for (p in seq_len(nperm)) {
    sc <- vapply(seq_len(n_contrasts), function(j) {
      n <- n_det_vec[j]
      if (is.na(n) || n == 0 || n > length(uni_vec[[j]])) return(NA_real_)
      mean(sample(uni_vec[[j]], n)) * log2(n + 1)
    }, numeric(1))
    out[p] <- mean(sc, na.rm = TRUE)
  }
  out
}

message("Running ", N_PERM, " permutations per axis and design...")
null_list <- vector("list", length(axes))
for (i in seq_along(axes)) {
  ax <- axes[i]
  K <- axis_info$n_axis_genes_in_universe[axis_info$regulator == ax]
  n_det_vec <- observed$n_detected[observed$regulator == ax][match(colnames(M), observed$contrast[observed$regulator == ax])]
  message("  axis: ", ax, " (K = ", K, ")")
  null_list[[i]] <- tibble(
    regulator = ax,
    perm_id = seq_len(N_PERM),
    null_score_common_set = null_common_set(K, N_PERM),
    null_score_per_contrast_draw = null_per_contrast(n_det_vec, N_PERM)
  )
}
null_df <- bind_rows(null_list)
readr::write_csv(null_df, file.path(table_dir, "regulator_axis_permutation_null_scores.csv"))

###############################################################################
# Summary with multiple-testing adjustment across the six axes
###############################################################################

emp_p <- function(null, obs) (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)

summary_df <- axis_info %>%
  rowwise() %>%
  mutate(
    null_mean_common_set = mean(null_df$null_score_common_set[null_df$regulator == regulator]),
    null_sd_common_set = sd(null_df$null_score_common_set[null_df$regulator == regulator]),
    empirical_p_common_set = emp_p(null_df$null_score_common_set[null_df$regulator == regulator], observed_mean_score),
    null_mean_per_contrast_draw = mean(null_df$null_score_per_contrast_draw[null_df$regulator == regulator]),
    null_sd_per_contrast_draw = sd(null_df$null_score_per_contrast_draw[null_df$regulator == regulator]),
    empirical_p_per_contrast_draw = emp_p(null_df$null_score_per_contrast_draw[null_df$regulator == regulator], observed_mean_score)
  ) %>%
  ungroup() %>%
  mutate(
    z_common_set = (observed_mean_score - null_mean_common_set) / null_sd_common_set,
    p_adj_BH_common_set = p.adjust(empirical_p_common_set, method = "BH"),
    p_adj_Holm_common_set = p.adjust(empirical_p_common_set, method = "holm"),
    p_adj_BH_per_contrast_draw = p.adjust(empirical_p_per_contrast_draw, method = "BH"),
    p_adj_Holm_per_contrast_draw = p.adjust(empirical_p_per_contrast_draw, method = "holm")
  ) %>%
  arrange(empirical_p_common_set)

readr::write_csv(summary_df, file.path(table_dir, "regulator_axis_permutation_summary.csv"))

message("Permutation test complete.")
cat("\n=== Regulatory-axis permutation test (", N_PERM, " permutations) ===\n", sep = "")
print(as.data.frame(summary_df %>% mutate(across(where(is.numeric), ~signif(.x, 4)))), row.names = FALSE)
cat("\nNF-kB observed mean score (authoritative):",
    round(summary_df$observed_mean_score[grepl("^NF-kB", summary_df$regulator)], 3), "\n")
