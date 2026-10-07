###############################################################################
# Tonsil complement/coagulation module (blue/ME2): robustness to lung-derived
# (alveolar/airway) transcripts
#
# Background (follows from 17_cell_composition_marker_scores.R):
#   Script 17 showed that SFTPC, SFTPB and AGER (alveolar markers with no
#   expected role in tonsil, and which are normally DOWNregulated, not
#   upregulated, in infected/injured lung) are near-absent in baseline tonsil
#   and uniformly high in every infected tonsil sample irrespective of DPI,
#   i.e. a baseline-vs-infected step change rather than a time trend. These
#   genes are themselves members of the tonsil blue/ME2 module that the
#   manuscript interprets as an infection-associated complement/coagulation
#   program, so the module-level association with infection may be partly
#   driven by lung-derived RNA in the infected tonsil samples.
#
#   This script asks whether the tonsil complement/coagulation signal survives
#   removal of that lung-derived signal, by:
#     (1) listing which alveolar/airway marker genes are module members and
#         their module membership (kME);
#     (2) recomputing the module eigengene after excluding ALL pre-specified
#         lung-specific marker genes found in the module, and comparing it to
#         the original ME2 (agreement, and association with infection/DPI);
#     (3) testing, gene by gene, whether key complement/coagulation genes
#         remain associated with infection after adjusting (partial Spearman
#         correlation) for a per-sample lung-contamination score.
#   The lung-gene list is fixed a priori (canonical type-II/type-I
#   pneumocyte and airway secretory markers), not selected from these data.
#
#   This is a robustness analysis. It cannot identify the SOURCE of the
#   lung-derived signal (sample handling, necropsy order, aspirated material,
#   batch): the public sample metadata contain no animal ID, batch, collection
#   method or necropsy order.
#
# Inputs:
#   advanced_analyses/tables/WGCNA_Tonsil_object.rds    (datExpr = VST, samples x genes; MEs)
#   advanced_analyses/tables/WGCNA_Tonsil_gene_modules.csv
#   GSE310471/results/tables/GSE310471_sample_metadata.csv
#
# Outputs:
#   advanced_analyses/tables/tonsil_ME2_lung_gene_membership.csv
#   advanced_analyses/tables/tonsil_ME2_original_vs_lung_excluded_summary.csv
#   advanced_analyses/tables/tonsil_ME2_per_sample_eigengenes.csv
#   advanced_analyses/tables/tonsil_complement_coagulation_genes_partial_correlation.csv
#   advanced_analyses/figures/tonsil_ME2_lung_gene_exclusion.png
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
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
gse_dir <- file.path(project_dir, "GSE310471", "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed a-priori lung-specific marker list (not data-derived).
#   Type-II pneumocyte: SFTPC, SFTPB, SFTPA1, SFTPA2, SFTPD, NAPSA, LAMP3,
#                       SLC34A2, ABCA3, LPCAT1
#   Type-I pneumocyte : AGER, HOPX, EMP2
#   Airway secretory  : SCGB1A1, SCGB3A2
#   Lung lineage TF   : NKX2-1
lung_genes <- c(
  "SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SFTPD", "NAPSA", "LAMP3", "SLC34A2",
  "ABCA3", "LPCAT1", "AGER", "HOPX", "EMP2", "SCGB1A1", "SCGB3A2", "NKX2-1"
)

# Core alveolar markers used to build the per-sample lung-contamination score
# (the unambiguous, high-specificity subset; HOPX/EMP2/LPCAT1/SLC34A2 are less
# lung-restricted and are excluded from the score but still excluded from the
# module when recomputing the eigengene).
contamination_score_genes <- c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "NAPSA", "AGER")

key_genes <- c(
  "C1QA", "C1QB", "C1QC", "C3", "CFB", "CFH", "CFI", "SERPING1", "C2", "C4A", "C4B",
  "FGA", "FGB", "FGG", "PROC", "PROS1", "THBD", "F3", "PLAU", "PLAUR", "SERPINE1", "VWF"
)

###############################################################################
# 1. Load tonsil WGCNA object, module assignment, sample metadata
###############################################################################

obj <- readRDS(file.path(table_dir, "WGCNA_Tonsil_object.rds"))
datExpr <- obj$datExpr          # samples x genes (VST)
MEs <- obj$MEs
gene_modules <- readr::read_csv(file.path(table_dir, "WGCNA_Tonsil_gene_modules.csv"), show_col_types = FALSE)

metadata <- readr::read_csv(file.path(gse_dir, "GSE310471_sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(
    dpi_numeric = ifelse(dpi == "baseline", 0, as.numeric(str_extract(dpi, "[0-9]+"))),
    infected = ifelse(dpi == "baseline", 0, 1)
  )

sample_meta <- tibble(sample_id = rownames(datExpr)) %>%
  left_join(metadata %>% select(sample_id, dpi, dpi_numeric, infected), by = "sample_id")
stopifnot(!anyNA(sample_meta$dpi_numeric))

module_num <- 2
me_name <- paste0("ME", module_num)
stopifnot(me_name %in% colnames(MEs))

module_genes <- gene_modules %>% filter(module_numeric == module_num) %>% pull(gene)
module_genes <- intersect(module_genes, colnames(datExpr))
cat("Module ", me_name, " (blue): ", length(module_genes), " genes in expression matrix\n", sep = "")

###############################################################################
# 2. Which lung-specific genes are in the module, and how strongly?
###############################################################################

kme_all <- as.numeric(cor(datExpr[, module_genes, drop = FALSE], MEs[[me_name]], use = "pairwise.complete.obs"))
names(kme_all) <- module_genes

lung_membership <- tibble(
  gene = lung_genes,
  detected_in_matrix = gene %in% colnames(datExpr),
  in_module_ME2 = gene %in% module_genes,
  assigned_module_numeric = gene_modules$module_numeric[match(gene, gene_modules$gene)],
  kME_vs_ME2 = ifelse(gene %in% module_genes, kme_all[gene], NA_real_)
)
lung_membership$kME_rank_in_module <- ifelse(
  lung_membership$in_module_ME2,
  rank(-abs(kme_all))[match(lung_membership$gene, names(kme_all))],
  NA_real_
)
readr::write_csv(lung_membership, file.path(table_dir, "tonsil_ME2_lung_gene_membership.csv"))

lung_in_module <- lung_membership %>% filter(in_module_ME2) %>% pull(gene)
cat("\nLung-specific genes in the module (", length(lung_in_module), " of ", length(module_genes), " module genes): ",
    paste(lung_in_module, collapse = ", "), "\n", sep = "")

###############################################################################
# 3. Recompute the module eigengene excluding lung-specific genes
###############################################################################

genes_kept <- setdiff(module_genes, lung_genes)
me_reduced <- moduleEigengenes(
  datExpr[, genes_kept, drop = FALSE],
  colors = rep("blue", length(genes_kept))
)$eigengenes[[1]]

# Orient sign to agree with the original ME2
if (cor(me_reduced, MEs[[me_name]]) < 0) me_reduced <- -me_reduced

# Re-run the original eigengene computation on the full module as a consistency check
me_full_recomputed <- moduleEigengenes(
  datExpr[, module_genes, drop = FALSE],
  colors = rep("blue", length(module_genes))
)$eigengenes[[1]]
if (cor(me_full_recomputed, MEs[[me_name]]) < 0) me_full_recomputed <- -me_full_recomputed

# Lung-contamination score: mean gene-wise z-score of core alveolar markers
core_present <- intersect(contamination_score_genes, colnames(datExpr))
z_core <- scale(datExpr[, core_present, drop = FALSE])
contamination_score <- rowMeans(z_core)

per_sample <- sample_meta %>%
  mutate(
    ME2_original = as.numeric(MEs[[me_name]]),
    ME2_full_recomputed = as.numeric(me_full_recomputed),
    ME2_lung_excluded = as.numeric(me_reduced),
    lung_contamination_score = as.numeric(contamination_score)
  )
readr::write_csv(per_sample, file.path(table_dir, "tonsil_ME2_per_sample_eigengenes.csv"))

spearman_row <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

summarise_me <- function(label, x) {
  a <- spearman_row(x, per_sample$dpi_numeric)
  b <- spearman_row(x, per_sample$infected)
  inf <- per_sample$infected == 1
  c3 <- spearman_row(x[inf], per_sample$dpi_numeric[inf])
  d <- spearman_row(x, per_sample$lung_contamination_score)
  tibble(
    eigengene = label,
    rho_vs_DPI = a["rho"], p_vs_DPI = a["p"],
    rho_vs_infected = b["rho"], p_vs_infected = b["p"],
    rho_vs_DPI_infected_only = c3["rho"], p_vs_DPI_infected_only = c3["p"],
    rho_vs_lung_contamination_score = d["rho"], p_vs_lung_contamination_score = d["p"]
  )
}

me_summary <- bind_rows(
  summarise_me("Original ME2 (all module genes)", per_sample$ME2_original),
  summarise_me("ME2 recomputed, lung genes excluded", per_sample$ME2_lung_excluded)
) %>%
  mutate(
    n_module_genes_used = c(length(module_genes), length(genes_kept)),
    rho_vs_original_ME2 = c(1, as.numeric(cor(per_sample$ME2_original, per_sample$ME2_lung_excluded, method = "spearman")))
  )
readr::write_csv(me_summary, file.path(table_dir, "tonsil_ME2_original_vs_lung_excluded_summary.csv"))

###############################################################################
# 4. Gene-level: do key complement/coagulation genes remain associated with
#    infection after adjusting for the lung-contamination score?
#    (partial Spearman: rank-transform, residualize on contamination rank)
###############################################################################

partial_spearman <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  ex <- resid(lm(rx ~ rz)); ey <- resid(lm(ry ~ rz))
  r <- cor(ex, ey)
  n <- length(x)
  tstat <- r * sqrt((n - 3) / (1 - r^2))
  p <- 2 * pt(-abs(tstat), df = n - 3)
  c(rho = r, p = p)
}

key_present <- intersect(key_genes, colnames(datExpr))
missing_key <- setdiff(key_genes, colnames(datExpr))

gene_results <- map_dfr(key_present, function(g) {
  x <- datExpr[, g]
  un <- spearman_row(x, per_sample$infected)
  adj <- partial_spearman(x, per_sample$infected, per_sample$lung_contamination_score)
  tibble(
    gene = g,
    module_numeric = gene_modules$module_numeric[match(g, gene_modules$gene)],
    rho_vs_infected_unadjusted = un["rho"], p_unadjusted = un["p"],
    rho_vs_infected_adj_lung_score = adj["rho"], p_adj_lung_score = adj["p"],
    rho_vs_lung_contamination_score = spearman_row(x, per_sample$lung_contamination_score)["rho"]
  )
}) %>%
  mutate(
    fdr_unadjusted = p.adjust(p_unadjusted, method = "BH"),
    fdr_adj_lung_score = p.adjust(p_adj_lung_score, method = "BH")
  ) %>%
  arrange(p_adj_lung_score)

readr::write_csv(gene_results, file.path(table_dir, "tonsil_complement_coagulation_genes_partial_correlation.csv"))

###############################################################################
# 5. Figure: original vs lung-excluded eigengene, and contamination vs DPI
###############################################################################

plot_long <- per_sample %>%
  mutate(dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI"))) %>%
  select(sample_id, dpi, ME2_original, ME2_lung_excluded, lung_contamination_score) %>%
  pivot_longer(
    c(ME2_original, ME2_lung_excluded, lung_contamination_score),
    names_to = "measure", values_to = "value"
  ) %>%
  mutate(measure = factor(
    measure,
    levels = c("ME2_original", "ME2_lung_excluded", "lung_contamination_score"),
    labels = c("Original ME2 (all module genes)", "ME2 recomputed without lung-specific genes",
               "Lung-contamination score (mean z of SFTPC/SFTPB/SFTPA1/SFTPA2/NAPSA/AGER)")
  ))

p <- ggplot(plot_long, aes(x = dpi, y = value)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
  facet_wrap(~ measure, ncol = 3, scales = "free_y",
             labeller = label_wrap_gen(width = 32)) +
  labs(
    title = "Tonsil blue/ME2 module with and without lung-specific genes",
    subtitle = "Eigengene recomputed after removing pre-specified alveolar/airway markers; lung-contamination score shown for reference",
    x = "Days post infection", y = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(figure_dir, "tonsil_ME2_lung_gene_exclusion.png"), p, width = 11, height = 4.2, dpi = 300)

###############################################################################
# 6. Console summary
###############################################################################

cat("\n=== Lung-specific genes: module membership ===\n")
print(as.data.frame(lung_membership))

cat("\n=== Module eigengene: original vs lung-genes-excluded ===\n")
print(as.data.frame(me_summary))
cat("\nConsistency check: Spearman rho between original ME2 and recomputed full-module eigengene: ",
    round(cor(per_sample$ME2_original, per_sample$ME2_full_recomputed, method = "spearman"), 3), "\n", sep = "")

cat("\n=== Key complement/coagulation genes: infection association, unadjusted vs adjusted for lung-contamination score ===\n")
print(as.data.frame(gene_results))
if (length(missing_key) > 0) cat("\nNot in the WGCNA input matrix (not among the 8000 most variable genes): ", paste(missing_key, collapse = ", "), "\n", sep = "")

n_sig_unadj <- sum(gene_results$fdr_unadjusted < 0.05)
n_sig_adj <- sum(gene_results$fdr_adj_lung_score < 0.05)
cat("\n", n_sig_unadj, " of ", nrow(gene_results), " key genes are FDR-significant for infection unadjusted; ",
    n_sig_adj, " remain FDR-significant after adjusting for the lung-contamination score.\n", sep = "")

cat("\nInterpretation guide:\n",
    "  - If ME2 recomputed without lung genes still tracks infection/DPI strongly (and agrees with the\n",
    "    original ME2), the complement/coagulation program is not explained by the alveolar genes.\n",
    "  - If the association collapses, the module signal in tonsil is largely attributable to lung-derived\n",
    "    transcripts and the tonsil complement/coagulation conclusion must be heavily qualified.\n",
    "  - Because the contamination score is itself perfectly confounded with infection in these data\n",
    "    (near zero at baseline, high in all infected samples), adjusting for it removes most infection\n",
    "    signal by construction; treat the adjusted gene-level result as a conservative lower bound and\n",
    "    prioritize the lung-genes-excluded eigengene result.\n", sep = "")

cat("\nOutputs written to advanced_analyses/tables and advanced_analyses/figures (tonsil_ME2_* and tonsil_complement_coagulation_* files).\n")
