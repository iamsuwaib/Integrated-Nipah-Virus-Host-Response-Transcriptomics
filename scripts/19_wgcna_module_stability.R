###############################################################################
# WGCNA disease-module stability analysis for GSE310471 lung and tonsil
#
# The main WGCNA analysis (03_wgcna_gse310471_exploratory.R) is exploratory by
# design given the modest per-tissue sample sizes (Lung n = 17, Tonsil n = 18).
# This script adds three complementary robustness checks for the two
# disease-associated modules the manuscript relies on (lung tan/ME12, tonsil
# blue/ME2), without altering the original module definitions, gene
# memberships, or any existing table/figure:
#
#   (A) Bootstrap stability of the module-trait correlation: resample samples
#       with replacement (same tissue, same frozen module gene list), recompute
#       the module eigengene and its correlation with infection status and
#       numeric DPI, and summarize the resulting distribution of r values.
#   (B) Leave-one-out (LOO) stability of the module-trait correlation: same
#       idea, but each of the n samples is dropped one at a time rather than
#       resampled.
#   (C) LOO network re-clustering: for each leave-one-out fold, the FULL
#       WGCNA pipeline (same soft power, same blockwiseModules parameters) is
#       re-run on all ~8000 genes with that one sample removed, and the
#       resulting module most similar to the original focal module (by gene-
#       set Jaccard overlap) is identified and checked for whether it remains
#       significantly correlated with the trait. This tests whether a similar
#       co-expression module re-emerges under sample perturbation, as a
#       lightweight resampling-based analog to formal module preservation
#       statistics (full permutation-based modulePreservation() Z-summary
#       statistics were not run given the computational cost and the already
#       modest per-tissue sample size).
#
# Inputs: WGCNA_Lung_object.rds and WGCNA_Tonsil_object.rds, saved by
# 03_wgcna_gse310471_exploratory.R (contain datExpr, net, MEs, traits,
# module_trait_cor, module_trait_p, soft_power for each tissue).
#
# Outputs (advanced_analyses/tables/):
#   WGCNA_stability_bootstrap.csv       - per-bootstrap r values, both traits, both tissues
#   WGCNA_stability_loo_eigengene.csv   - per-LOO-fold r values (frozen module), both traits
#   WGCNA_stability_loo_reclustering.csv- per-LOO-fold Jaccard overlap + best-match r
#   WGCNA_stability_summary.csv         - one-row-per-tissue summary of all three checks
#
# Runtime note: parts (A) and (B) are fast (seconds). Part (C) re-runs
# blockwiseModules once per left-out sample (17 + 18 = 35 total reruns on the
# full ~8000-gene expression matrix) and is the slow part of this script -
# expect roughly 20-60 minutes total depending on the machine, similar in
# per-run cost to the original WGCNA run.
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "WGCNA")
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

suppressPackageStartupMessages({
  library(tidyverse)
  library(WGCNA)
})

allowWGCNAThreads()

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")

set.seed(20260829)
N_BOOT <- 1000

lung_obj <- readRDS(file.path(table_dir, "WGCNA_Lung_object.rds"))
tonsil_obj <- readRDS(file.path(table_dir, "WGCNA_Tonsil_object.rds"))

frozen_module_eigengene_cor <- function(sub_expr, sub_trait) {
  # sub_expr: samples x genes matrix for the frozen focal-module gene set only
  rownames(sub_expr) <- paste0("s", seq_len(nrow(sub_expr)))
  me <- tryCatch(
    as.numeric(moduleEigengenes(sub_expr, colors = rep("focal", ncol(sub_expr)))$eigengenes[, 1]),
    error = function(e) rep(NA_real_, nrow(sub_expr))
  )
  # as.numeric(...)[1] guards against cor() ever returning a 1x1 matrix
  # instead of a plain scalar, which otherwise turns the result column into
  # an unwritable list-column when thousands of rows are combined.
  as.numeric(suppressWarnings(cor(me, as.numeric(sub_trait), use = "pairwise.complete.obs")))[1]
}

run_stability <- function(obj, focal_module_num, focal_module_label, traits_to_test = c("infected", "dpi_numeric")) {
  tissue <- obj$tissue
  datExpr <- obj$datExpr
  net <- obj$net
  traits <- obj$traits
  soft_power <- obj$soft_power
  n <- nrow(datExpr)

  focal_genes <- colnames(datExpr)[net$colors == focal_module_num]
  message(tissue, " focal module ", focal_module_label, ": ", length(focal_genes), " genes")

  ## --- (A) Bootstrap stability of frozen-module eigengene trait correlation ---
  boot_rows <- vector("list", N_BOOT * length(traits_to_test))
  k <- 1
  for (b in seq_len(N_BOOT)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    sub_expr <- datExpr[idx, focal_genes, drop = FALSE]
    for (tr in traits_to_test) {
      sub_trait <- traits[[tr]][idx]
      r <- frozen_module_eigengene_cor(sub_expr, sub_trait)
      boot_rows[[k]] <- tibble(tissue = tissue, module = focal_module_label, trait = tr, boot_id = b, r = r)
      k <- k + 1
    }
  }
  boot_df <- bind_rows(boot_rows)

  ## --- (B) LOO stability of frozen-module eigengene trait correlation ---
  loo_rows <- vector("list", n * length(traits_to_test))
  k <- 1
  for (i in seq_len(n)) {
    sub_expr <- datExpr[-i, focal_genes, drop = FALSE]
    for (tr in traits_to_test) {
      sub_trait <- traits[[tr]][-i]
      r <- frozen_module_eigengene_cor(sub_expr, sub_trait)
      loo_rows[[k]] <- tibble(
        tissue = tissue, module = focal_module_label, trait = tr,
        left_out_index = i, left_out_sample = rownames(datExpr)[i], r = r
      )
      k <- k + 1
    }
  }
  loo_df <- bind_rows(loo_rows)

  ## --- (C) LOO full network re-clustering ---
  reclust_rows <- vector("list", n)
  for (i in seq_len(n)) {
    message(tissue, " LOO re-clustering fold ", i, " of ", n)
    sub_expr_full <- datExpr[-i, , drop = FALSE]
    sub_net <- blockwiseModules(
      sub_expr_full,
      power = soft_power,
      networkType = "signed",
      TOMType = "signed",
      minModuleSize = 30,
      reassignThreshold = 0,
      mergeCutHeight = 0.25,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    )
    genes_present <- colnames(sub_expr_full)
    mod_ids <- sort(unique(sub_net$colors))
    overlaps <- sapply(mod_ids, function(m) {
      mod_genes <- genes_present[sub_net$colors == m]
      length(intersect(mod_genes, focal_genes)) / length(union(mod_genes, focal_genes))
    })
    best_idx <- which.max(overlaps)
    best_mod_id <- mod_ids[best_idx]
    best_jaccard <- overlaps[best_idx]

    MEs_sub <- tryCatch(
      moduleEigengenes(sub_expr_full, colors = sub_net$colors)$eigengenes,
      error = function(e) NULL
    )
    r_infected <- NA_real_
    r_dpi <- NA_real_
    if (!is.null(MEs_sub)) {
      me_name <- paste0("ME", best_mod_id)
      if (me_name %in% colnames(MEs_sub)) {
        r_infected <- as.numeric(suppressWarnings(cor(as.numeric(MEs_sub[[me_name]]), as.numeric(traits$infected[-i]), use = "p")))[1]
        r_dpi <- as.numeric(suppressWarnings(cor(as.numeric(MEs_sub[[me_name]]), as.numeric(traits$dpi_numeric[-i]), use = "p")))[1]
      }
    }
    reclust_rows[[i]] <- tibble(
      tissue = tissue, module = focal_module_label,
      left_out_index = i, left_out_sample = rownames(datExpr)[i],
      best_match_jaccard = best_jaccard,
      n_genes_best_match = sum(sub_net$colors == best_mod_id),
      r_infected_best_match = r_infected,
      r_dpi_numeric_best_match = r_dpi
    )
  }
  reclust_df <- bind_rows(reclust_rows)

  list(boot = boot_df, loo = loo_df, reclust = reclust_df)
}

lung_stab <- run_stability(lung_obj, focal_module_num = 12, focal_module_label = "tan/ME12")
tonsil_stab <- run_stability(tonsil_obj, focal_module_num = 2, focal_module_label = "blue/ME2")

boot_all <- bind_rows(lung_stab$boot, tonsil_stab$boot)
loo_all <- bind_rows(lung_stab$loo, tonsil_stab$loo)
reclust_all <- bind_rows(lung_stab$reclust, tonsil_stab$reclust)

readr::write_csv(boot_all, file.path(table_dir, "WGCNA_stability_bootstrap.csv"))
readr::write_csv(loo_all, file.path(table_dir, "WGCNA_stability_loo_eigengene.csv"))
readr::write_csv(reclust_all, file.path(table_dir, "WGCNA_stability_loo_reclustering.csv"))

summary_df <- bind_rows(
  boot_all %>%
    group_by(tissue, module, trait) %>%
    summarise(
      boot_mean_r = mean(r, na.rm = TRUE),
      boot_sd_r = sd(r, na.rm = TRUE),
      boot_q2.5 = quantile(r, 0.025, na.rm = TRUE),
      boot_q97.5 = quantile(r, 0.975, na.rm = TRUE),
      boot_pct_same_sign_as_original = mean(sign(r) == sign(boot_mean_r), na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(-c(tissue, module, trait), names_to = "metric", values_to = "value"),
  loo_all %>%
    group_by(tissue, module, trait) %>%
    summarise(
      loo_mean_r = mean(r, na.rm = TRUE),
      loo_min_r = min(r, na.rm = TRUE),
      loo_max_r = max(r, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-c(tissue, module, trait), names_to = "metric", values_to = "value")
)
readr::write_csv(summary_df, file.path(table_dir, "WGCNA_stability_summary.csv"))

reclust_summary <- reclust_all %>%
  group_by(tissue, module) %>%
  summarise(
    mean_jaccard = mean(best_match_jaccard, na.rm = TRUE),
    min_jaccard = min(best_match_jaccard, na.rm = TRUE),
    mean_r_infected = mean(r_infected_best_match, na.rm = TRUE),
    mean_r_dpi = mean(r_dpi_numeric_best_match, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(reclust_summary, file.path(table_dir, "WGCNA_stability_reclustering_summary.csv"))

message("WGCNA module stability analysis complete.")
print(summary_df)
print(reclust_summary)
