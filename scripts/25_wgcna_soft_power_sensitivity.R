###############################################################################
# WGCNA soft-thresholding power: selection rule and sensitivity analysis
#
# Purpose
#   (1) Document the soft-power selection rule used in
#       07_wgcna_gse310471_exploratory.R. The rule is:
#         soft_power <- pickSoftThreshold(..., networkType = "signed")$powerEstimate
#         (lowest tested power whose signed scale-free fit R^2 reaches the
#          pickSoftThreshold default RsquaredCut = 0.85);
#         if powerEstimate is NA (no tested power reaches the cut), a fixed
#         pre-specified fallback of 6 is used.
#       Lung reached the criterion (power 18); tonsil did not (maximum signed
#       fit ~0.55 at power 20), so the tonsil network used the fallback of 6.
#   (2) Test whether the key biological interpretation depends on that choice.
#       The identical network construction used in script 07 (same expression
#       matrix, signed network, signed TOM, minModuleSize = 30, mergeCutHeight
#       = 0.25, reassignThreshold = 0) is re-run over a range of powers, and for
#       each power the module that best matches the original focal module
#       (lung tan/ME12; tonsil blue/ME2) is identified by gene-set Jaccard
#       overlap. For each power we report:
#         - scale-free fit (signed R^2) and mean connectivity,
#         - size and Jaccard overlap of the best-matching module,
#         - the fraction of the original module's genes retained,
#         - how many of the displayed core candidate genes (Figure 4A) fall in
#           the best-matching module and in the single module holding most of them,
#         - eigengene correlation of the best-matching module with infection
#           status and numeric DPI.
#       Infection status and numeric DPI are strongly related traits in this
#       design (baseline = 0, all infected samples > 0); they are reported side by
#       side for completeness but are not independent confirmations.
#
# Inputs (written by earlier scripts):
#   advanced_analyses/tables/WGCNA_Lung_object.rds, WGCNA_Tonsil_object.rds
#   advanced_analyses/tables/focused_wgcna_candidate_modules.csv
# Outputs (advanced_analyses/tables, advanced_analyses/figures):
#   WGCNA_soft_power_sensitivity.csv
#   WGCNA_soft_power_scale_free_fit.csv
#   WGCNA_soft_power_sensitivity.png
# Runtime: one blockwiseModules run per tissue-power combination (about
# 17 runs on ~8000 genes); expect roughly 20-40 minutes.
###############################################################################

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(WGCNA)
})
allowWGCNAThreads()

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")

lung_obj <- readRDS(file.path(table_dir, "WGCNA_Lung_object.rds"))
tonsil_obj <- readRDS(file.path(table_dir, "WGCNA_Tonsil_object.rds"))
focused <- readr::read_csv(file.path(table_dir, "focused_wgcna_candidate_modules.csv"), show_col_types = FALSE)

# sanity check of module identity (numeric label -> colour) in the original runs
stopifnot(labels2colors(12) == "tan", labels2colors(2) == "blue")

powers_to_test <- list(
  Lung   = c(6, 12, 14, 16, 18, 20),
  Tonsil = c(4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20)
)
focal_label <- c(Lung = "tan/ME12", Tonsil = "blue/ME2")
focal_num   <- c(Lung = 12, Tonsil = 2)

fit_rows <- list()
res_rows <- list()

for (obj in list(lung_obj, tonsil_obj)) {
  tissue <- obj$tissue
  datExpr <- obj$datExpr
  traits <- obj$traits
  focal_genes <- colnames(datExpr)[obj$net$colors == focal_num[[tissue]]]
  tis_name <- tissue
  focal_n <- focal_num[[tis_name]]
  core_genes <- focused %>%
    filter(.data$tissue == tis_name, .data$module_numeric == focal_n) %>%
    pull(gene) %>% unique()
  message(tissue, ": original focal module ", length(focal_genes), " genes; core candidates ", length(core_genes))

  # scale-free fit over the same grid as script 07 plus the extra powers tested here
  grid <- sort(unique(c(1:10, seq(12, 20, by = 2), powers_to_test[[tissue]])))
  sft <- pickSoftThreshold(datExpr, powerVector = grid, verbose = 0, networkType = "signed")
  fit_rows[[tissue]] <- sft$fitIndices %>%
    transmute(tissue = tis_name, power = Power,
              signed_R2 = -sign(slope) * SFT.R.sq, slope = slope,
              mean_connectivity = sft$fitIndices[, 5],
              powerEstimate_default_0.85 = sft$powerEstimate)

  for (pw in powers_to_test[[tissue]]) {
    message(tissue, " power ", pw)
    net <- blockwiseModules(
      datExpr, power = pw, networkType = "signed", TOMType = "signed",
      minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
      numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0
    )
    genes_all <- colnames(datExpr)
    mods <- sort(setdiff(unique(net$colors), 0))
    jac <- sapply(mods, function(m) {
      g <- genes_all[net$colors == m]
      length(intersect(g, focal_genes)) / length(union(g, focal_genes))
    })
    best <- mods[which.max(jac)]
    best_genes <- genes_all[net$colors == best]
    # module that holds most of the core candidates (may differ from best match)
    core_mod_tab <- table(net$colors[match(core_genes, genes_all)])
    core_mod_tab <- core_mod_tab[names(core_mod_tab) != "0"]
    modal_mod <- if (length(core_mod_tab)) as.integer(names(core_mod_tab)[which.max(core_mod_tab)]) else NA_integer_
    modal_n <- if (length(core_mod_tab)) max(core_mod_tab) else 0L

    MEs <- moduleEigengenes(datExpr, colors = net$colors)$eigengenes
    me <- MEs[[paste0("ME", best)]]
    ct_inf <- suppressWarnings(cor.test(me, traits$infected))
    ct_dpi <- suppressWarnings(cor.test(me, traits$dpi_numeric))

    res_rows[[length(res_rows) + 1]] <- tibble(
      tissue = tis_name, module = focal_label[[tis_name]], power = pw,
      is_original_power = pw == obj$soft_power,
      n_modules = length(mods), n_unassigned_grey = sum(net$colors == 0),
      best_match_module = best, best_match_n_genes = length(best_genes),
      jaccard_with_original = max(jac),
      frac_original_genes_retained = length(intersect(best_genes, focal_genes)) / length(focal_genes),
      n_core_candidates = length(core_genes),
      n_core_in_best_match = sum(core_genes %in% best_genes),
      n_core_in_single_modal_module = modal_n,
      modal_module_equals_best_match = identical(as.integer(modal_mod), as.integer(best)),
      eigengene_r_infected = unname(ct_inf$estimate), eigengene_p_infected = ct_inf$p.value,
      eigengene_r_dpi_numeric = unname(ct_dpi$estimate), eigengene_p_dpi_numeric = ct_dpi$p.value
    )
  }
}

fit_df <- bind_rows(fit_rows)
res_df <- bind_rows(res_rows)
readr::write_csv(fit_df, file.path(table_dir, "WGCNA_soft_power_scale_free_fit.csv"))
readr::write_csv(res_df, file.path(table_dir, "WGCNA_soft_power_sensitivity.csv"))

# check: at the original power the re-run should reproduce the original module
chk <- res_df %>% filter(is_original_power) %>% select(tissue, power, jaccard_with_original)
print(chk)

p <- res_df %>%
  mutate(label = paste(tissue, module)) %>%
  pivot_longer(c(jaccard_with_original, eigengene_r_infected),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(power, value, colour = metric)) +
  geom_line() + geom_point() +
  facet_wrap(~ label, scales = "free_x") +
  labs(x = "Soft-thresholding power", y = "Value", colour = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(figure_dir, "WGCNA_soft_power_sensitivity.png"), p, width = 9, height = 4, dpi = 300)

message("Done")
print(res_df %>% select(tissue, power, best_match_n_genes, jaccard_with_original,
                        frac_original_genes_retained, n_core_in_best_match, n_core_candidates,
                        eigengene_r_infected, eigengene_r_dpi_numeric), n = Inf)
