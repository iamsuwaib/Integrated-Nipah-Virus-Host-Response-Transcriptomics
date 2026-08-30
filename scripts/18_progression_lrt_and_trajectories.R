###############################################################################
# GSE310471 temporal-progression analysis: formal DESeq2 LRT + module
# eigengene trajectories across DPI
#
# Purpose:
#   "Progression" was previously supported only by increasing DEG counts
#   across 3/4/5 DPI versus baseline, which is not a formal time/trend test.
#   A more rigorous characterization requires either (i) a formal time/trend
#   statistical model, or (ii) gene/module trajectories with uncertainty
#   estimates.
#
#   This script provides both:
#     (1) A formal DESeq2 likelihood-ratio test (LRT) per tissue, comparing
#         the full model (~dpi, i.e. the same per-tissue model already used
#         for the pairwise DPI-vs-baseline Wald contrasts) against a reduced
#         model with no dpi term (~1). This is a single omnibus test per gene
#         for whether expression varies at all across the four DPI levels,
#         independent of any specific pairwise contrast or significance
#         threshold used elsewhere in this study.
#     (2) Trajectory plots of the lung antiviral (tan/ME12) and tonsil
#         complement/coagulation (blue/ME2) WGCNA module eigengenes across
#         DPI, showing per-sample dispersion (the "uncertainty estimate").
#         The module-eigengene-vs-DPI correlations themselves are already
#         reported in the manuscript from the original WGCNA module-trait
#         analysis (Methods 2.5) and are not recomputed here.
#
# Inputs:
#   GSE310471/results/GSE310471_analysis_objects.rds (per-tissue fitted dds)
#   GSE310471/results/tables/GSE310471_sample_metadata.csv
#   advanced_analyses/tables/WGCNA_Lung_object.rds   (for ME12)
#   advanced_analyses/tables/WGCNA_Tonsil_object.rds (for ME2)
#
# Outputs:
#   advanced_analyses/tables/GSE310471_progression_LRT_full.csv
#   advanced_analyses/tables/GSE310471_progression_LRT_summary.csv
#   advanced_analyses/figures/progression_module_eigengene_trajectories.png
#   Console summary
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "DESeq2", "ggplot2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "DESeq2") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
})

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
gse_dir <- file.path(project_dir, "GSE310471")
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. Formal DESeq2 likelihood-ratio test (LRT) per tissue
###############################################################################

analysis_objects <- readRDS(file.path(gse_dir, "results", "GSE310471_analysis_objects.rds"))

run_lrt <- function(tissue_name) {
  message("Running DESeq2 LRT for tissue: ", tissue_name)
  dds <- analysis_objects$tissue_results[[tissue_name]]$dds
  # Same fitted per-tissue dds object used for the pairwise Wald contrasts
  # reported elsewhere (Methods 2.2, design = ~dpi); re-fit here with the
  # LRT test comparing the full model against a dpi-free reduced model.
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1, quiet = TRUE)
  res <- results(dds_lrt) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(tissue = tissue_name)
  res
}

lrt_lung <- run_lrt("Lung")
lrt_tonsil <- run_lrt("Tonsil")
lrt_all <- dplyr::bind_rows(lrt_lung, lrt_tonsil)

readr::write_csv(lrt_all, file.path(table_dir, "GSE310471_progression_LRT_full.csv"))

lrt_summary <- lrt_all %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(
    n_genes_tested = dplyr::n(),
    n_sig_LRT_padj_0.05 = sum(padj < 0.05, na.rm = TRUE),
    pct_sig_LRT = round(100 * n_sig_LRT_padj_0.05 / n_genes_tested, 1),
    .groups = "drop"
  )

readr::write_csv(lrt_summary, file.path(table_dir, "GSE310471_progression_LRT_summary.csv"))

###############################################################################
# 2. Module eigengene trajectories across DPI (lung tan/ME12, tonsil blue/ME2)
###############################################################################

metadata <- readr::read_csv(file.path(gse_dir, "results", "tables", "GSE310471_sample_metadata.csv"), show_col_types = FALSE)

load_module_eigengene <- function(rds_path, me_col) {
  obj <- readRDS(rds_path)
  MEs <- obj$MEs
  if (!me_col %in% colnames(MEs)) {
    stop("Module eigengene column '", me_col, "' not found in ", rds_path,
         ". Available columns: ", paste(colnames(MEs), collapse = ", "))
  }
  tibble::tibble(sample_id = rownames(MEs), module_eigengene = MEs[[me_col]])
}

lung_me <- load_module_eigengene(file.path(table_dir, "WGCNA_Lung_object.rds"), "ME12") %>%
  dplyr::left_join(metadata, by = "sample_id") %>%
  dplyr::mutate(tissue = "Lung", module_label = "tan/ME12 (lung antiviral module)")

tonsil_me <- load_module_eigengene(file.path(table_dir, "WGCNA_Tonsil_object.rds"), "ME2") %>%
  dplyr::left_join(metadata, by = "sample_id") %>%
  dplyr::mutate(tissue = "Tonsil", module_label = "blue/ME2 (tonsil complement/coagulation module)")

traj_df <- dplyr::bind_rows(lung_me, tonsil_me) %>%
  dplyr::mutate(dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI")))

readr::write_csv(
  traj_df %>% dplyr::select(sample_id, tissue, module_label, dpi, module_eigengene),
  file.path(table_dir, "GSE310471_progression_module_eigengene_trajectories.csv")
)

p <- ggplot(traj_df, aes(x = dpi, y = module_eigengene, color = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
  facet_wrap(~module_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Disease-module eigengene trajectories across infection time",
    subtitle = "Lung tan/ME12 (antiviral) and tonsil blue/ME2 (complement/coagulation)\nmodule eigengenes by days post infection",
    x = "Days post infection",
    y = "Module eigengene"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "progression_module_eigengene_trajectories.png"),
  p, width = 7.5, height = 8.5, dpi = 300
)

###############################################################################
# Console summary
###############################################################################

cat("\n=== Progression LRT summary ===\n")
print(as.data.frame(lrt_summary))
cat("\nFull LRT results written to: GSE310471_progression_LRT_full.csv\n")
cat("Trajectory data written to: GSE310471_progression_module_eigengene_trajectories.csv\n")
cat("Trajectory plot written to: progression_module_eigengene_trajectories.png\n")
