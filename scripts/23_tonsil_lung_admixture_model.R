###############################################################################
# Tonsil lung-RNA admixture model: can lung-derived RNA explain the infection-
# associated complement/coagulation signal in tonsil?
#
# Background (follows from scripts 17 and 22):
#   Script 22 showed that (i) only 5 of 1490 blue/ME2 module genes are lung-
#   specific, so removing them does not change the module eigengene, and (ii)
#   adjusting individual genes for a lung-contamination score removes their
#   infection association, but that adjustment is not diagnostic here because
#   the score is almost perfectly confounded with infection status. Both
#   results therefore leave the key question open: is the complement/
#   coagulation increase in infected tonsil more than lung-derived RNA could
#   plausibly contribute?
#
#   This script answers that with a simple, explicit mixing model instead of
#   a correlation adjustment:
#     1. A lung reference expression profile is taken from the GSE310471 lung
#        samples (mean normalized count per gene).
#     2. For each tonsil sample, the fraction f of lung-derived RNA is
#        estimated from lung-specific markers as the median across markers of
#        (tonsil normalized count / lung reference normalized count). Lung-
#        specific markers are essentially absent from uncontaminated tonsil,
#        so this ratio approximates the admixture fraction. The spread of the
#        per-marker estimates is reported as a check that a single fraction
#        describes the data.
#     3. The count of any gene G expected from this admixture alone is
#        f * (lung reference count of G). Subtracting it from the observed
#        tonsil count gives an admixture-corrected count.
#     4. For each gene, the observed infected-vs-baseline increase is compared
#        with the expected admixture contribution (explained fraction), and
#        the infection association is re-tested on admixture-corrected counts.
#
#   Interpretation: genes whose observed increase greatly exceeds the
#   expected admixture contribution remain supported as tonsil-intrinsic
#   changes. Genes whose increase is comparable to, or smaller than, the
#   expected contribution cannot be distinguished from lung-derived RNA.
#
# Assumptions and limits (to be stated in the manuscript):
#   - Tonsil and lung counts were normalized separately (per-tissue DESeq2
#     size factors), so f is an approximate, not an absolute, RNA fraction;
#     results are interpreted as orders of magnitude.
#   - Infection changes lung expression, so lung reference choice matters;
#     both an all-lung and a baseline-lung reference are run.
#   - Admixture of other tissues (e.g. blood, airway) is not modelled; the
#     public metadata contain no animal, batch, collection or necropsy
#     information, so the SOURCE of the lung signal cannot be determined.
#
# Inputs:
#   GSE310471/results/tables/GSE310471_Lung_normalized_counts.csv
#   GSE310471/results/tables/GSE310471_Tonsil_normalized_counts.csv
#   GSE310471/results/tables/GSE310471_sample_metadata.csv
#   advanced_analyses/tables/WGCNA_Tonsil_gene_modules.csv
#
# Outputs:
#   advanced_analyses/tables/tonsil_lung_admixture_per_sample_fraction.csv
#   advanced_analyses/tables/tonsil_lung_admixture_per_marker_fraction.csv
#   advanced_analyses/tables/tonsil_lung_admixture_gene_results.csv
#   advanced_analyses/tables/tonsil_lung_admixture_summary.csv
#   advanced_analyses/figures/tonsil_lung_admixture_explained_fraction.png
#   Console summary
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
gse_dir <- file.path(project_dir, "GSE310471", "results", "tables")
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

marker_genes <- c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SFTPD", "NAPSA", "AGER", "NKX2-1", "LAMP3", "ABCA3")
min_lung_reference <- 100   # lung reference normalized count required for a marker to be informative

key_genes <- c(
  "C1QA", "C1QB", "C1QC", "C3", "CFB", "CFH", "CFI", "SERPING1", "C2", "C4A", "C4B",
  "FGA", "FGB", "FGG", "PROC", "PROS1", "THBD", "F3", "PLAU", "PLAUR", "SERPINE1", "VWF"
)

###############################################################################
# 1. Load data
###############################################################################

read_counts <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  m <- x %>% column_to_rownames("Gene") %>% as.matrix()
  m
}

lung <- read_counts(file.path(gse_dir, "GSE310471_Lung_normalized_counts.csv"))
tonsil <- read_counts(file.path(gse_dir, "GSE310471_Tonsil_normalized_counts.csv"))

metadata <- readr::read_csv(file.path(gse_dir, "GSE310471_sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(infected = ifelse(dpi == "baseline", 0, 1))

tonsil_meta <- tibble(sample_id = colnames(tonsil)) %>%
  left_join(metadata %>% select(sample_id, dpi, infected), by = "sample_id")
lung_meta <- tibble(sample_id = colnames(lung)) %>%
  left_join(metadata %>% select(sample_id, dpi, infected), by = "sample_id")
stopifnot(!anyNA(tonsil_meta$infected), !anyNA(lung_meta$infected))

common_genes <- intersect(rownames(lung), rownames(tonsil))
cat("Genes shared between lung and tonsil matrices: ", length(common_genes), "\n", sep = "")

gene_modules <- readr::read_csv(file.path(table_dir, "WGCNA_Tonsil_gene_modules.csv"), show_col_types = FALSE)

###############################################################################
# 2. Fit the admixture model for each lung-reference choice
###############################################################################

spearman_p <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

run_model <- function(ref_label, ref_samples) {
  lung_ref <- rowMeans(lung[common_genes, ref_samples, drop = FALSE])
  tons <- tonsil[common_genes, , drop = FALSE]

  markers <- intersect(marker_genes, common_genes)
  markers <- markers[lung_ref[markers] >= min_lung_reference]
  stopifnot(length(markers) >= 3)

  # per-marker, per-sample admixture fraction estimates
  ratio_mat <- sweep(tons[markers, , drop = FALSE], 1, lung_ref[markers], "/")
  f_sample <- apply(ratio_mat, 2, median)

  per_marker <- as_tibble(ratio_mat, rownames = "marker") %>%
    pivot_longer(-marker, names_to = "sample_id", values_to = "fraction") %>%
    left_join(tonsil_meta, by = "sample_id") %>%
    mutate(reference = ref_label)

  per_sample <- tonsil_meta %>%
    mutate(lung_fraction_estimate = as.numeric(f_sample[sample_id]), reference = ref_label)

  # admixture-corrected counts
  expected <- outer(lung_ref, f_sample)                       # genes x samples
  corrected <- pmax(tons - expected[rownames(tons), colnames(tons)], 0)

  inf <- tonsil_meta$infected == 1
  base_mean <- rowMeans(tons[, !inf, drop = FALSE])
  inf_mean <- rowMeans(tons[, inf, drop = FALSE])
  exp_inf_mean <- rowMeans(expected[rownames(tons), inf, drop = FALSE])

  gene_res <- map_dfr(rownames(tons), function(g) {
    raw <- spearman_p(log2(tons[g, ] + 1), tonsil_meta$infected)
    cor_res <- spearman_p(log2(corrected[g, ] + 1), tonsil_meta$infected)
    tibble(
      gene = g,
      tonsil_baseline_mean = base_mean[g],
      tonsil_infected_mean = inf_mean[g],
      lung_reference = lung_ref[g],
      observed_increase = inf_mean[g] - base_mean[g],
      expected_admixture_contribution = exp_inf_mean[g],
      rho_infected_raw = raw["rho"], p_raw = raw["p"],
      rho_infected_corrected = cor_res["rho"], p_corrected = cor_res["p"]
    )
  }) %>%
    mutate(
      explained_fraction = ifelse(observed_increase > 0, expected_admixture_contribution / observed_increase, NA_real_),
      fdr_raw = p.adjust(p_raw, method = "BH"),
      fdr_corrected = p.adjust(p_corrected, method = "BH"),
      reference = ref_label,
      module_numeric = gene_modules$module_numeric[match(gene, gene_modules$gene)]
    )

  list(per_sample = per_sample, per_marker = per_marker, gene_res = gene_res, markers = markers)
}

lung_all <- colnames(lung)
lung_base <- lung_meta$sample_id[lung_meta$infected == 0]

res_all <- run_model("lung_all_samples", lung_all)
res_base <- run_model("lung_baseline_only", lung_base)

per_sample_out <- bind_rows(res_all$per_sample, res_base$per_sample)
per_marker_out <- bind_rows(res_all$per_marker, res_base$per_marker)
gene_out <- bind_rows(res_all$gene_res, res_base$gene_res)

readr::write_csv(per_sample_out, file.path(table_dir, "tonsil_lung_admixture_per_sample_fraction.csv"))
readr::write_csv(per_marker_out, file.path(table_dir, "tonsil_lung_admixture_per_marker_fraction.csv"))
readr::write_csv(gene_out, file.path(table_dir, "tonsil_lung_admixture_gene_results.csv"))

###############################################################################
# 3. Summaries: key genes and the blue/ME2 module as a whole
###############################################################################

summarise_reference <- function(res, ref_label) {
  g <- res$gene_res
  mod <- g %>% filter(module_numeric == 2)
  up_raw <- mod %>% filter(fdr_raw < 0.05, observed_increase > 0)
  still <- up_raw %>% filter(fdr_corrected < 0.05, rho_infected_corrected > 0)
  tibble(
    reference = ref_label,
    markers_used = paste(res$markers, collapse = ", "),
    median_fraction_baseline = median(res$per_sample$lung_fraction_estimate[res$per_sample$infected == 0]),
    median_fraction_infected = median(res$per_sample$lung_fraction_estimate[res$per_sample$infected == 1]),
    min_fraction_infected = min(res$per_sample$lung_fraction_estimate[res$per_sample$infected == 1]),
    max_fraction_infected = max(res$per_sample$lung_fraction_estimate[res$per_sample$infected == 1]),
    n_module_genes_tested = nrow(mod),
    n_module_genes_up_in_infected_raw = nrow(up_raw),
    n_of_those_still_up_after_admixture_correction = nrow(still),
    median_explained_fraction_up_module_genes = median(up_raw$explained_fraction, na.rm = TRUE),
    n_up_module_genes_with_explained_fraction_gt_0.5 = sum(up_raw$explained_fraction > 0.5, na.rm = TRUE)
  )
}

summary_out <- bind_rows(
  summarise_reference(res_all, "lung_all_samples"),
  summarise_reference(res_base, "lung_baseline_only")
)
readr::write_csv(summary_out, file.path(table_dir, "tonsil_lung_admixture_summary.csv"))

key_out <- gene_out %>%
  filter(gene %in% key_genes) %>%
  select(reference, gene, module_numeric, tonsil_baseline_mean, tonsil_infected_mean, lung_reference,
         observed_increase, expected_admixture_contribution, explained_fraction,
         rho_infected_raw, fdr_raw, rho_infected_corrected, fdr_corrected) %>%
  arrange(reference, desc(observed_increase))

###############################################################################
# 4. Figure: explained fraction for module genes up in infected tonsil
###############################################################################

plot_df <- gene_out %>%
  filter(reference == "lung_all_samples", module_numeric == 2, fdr_raw < 0.05, observed_increase > 0) %>%
  mutate(
    label = ifelse(gene %in% key_genes, gene, NA_character_),
    explained_plot = pmin(explained_fraction, 2)
  )

p <- ggplot(plot_df, aes(x = log10(observed_increase + 1), y = explained_plot)) +
  geom_hline(yintercept = c(0.5, 1), linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_point(data = plot_df %>% filter(!is.na(label)), colour = "firebrick", size = 2.2) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 40, na.rm = TRUE) +
  labs(
    title = "Tonsil blue/ME2 genes up in infected tonsil: fraction of the increase expected from lung-RNA admixture",
    subtitle = "Dashed lines: 50% and 100% of the observed increase explained by the modelled admixture (values capped at 2)",
    x = "log10(observed increase in normalized count + 1)",
    y = "Expected admixture contribution / observed increase"
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(figure_dir, "tonsil_lung_admixture_explained_fraction.png"), p, width = 9, height = 6.5, dpi = 300)

###############################################################################
# 5. Console summary
###############################################################################

cat("\n=== Estimated lung-derived RNA fraction in tonsil (per sample) ===\n")
print(as.data.frame(per_sample_out %>% filter(reference == "lung_all_samples") %>%
                      arrange(infected, dpi) %>% select(sample_id, dpi, lung_fraction_estimate)))

cat("\n=== Per-marker agreement of the fraction estimate (infected samples, lung_all_samples reference) ===\n")
print(as.data.frame(
  per_marker_out %>% filter(reference == "lung_all_samples", infected == 1) %>%
    group_by(marker) %>% summarise(median_fraction = median(fraction), min = min(fraction), max = max(fraction), .groups = "drop")
))

cat("\n=== Model summary ===\n")
print(as.data.frame(summary_out))

cat("\n=== Key complement/coagulation genes (observed increase vs expected admixture contribution) ===\n")
print(as.data.frame(key_out), digits = 3)

cat("\nInterpretation guide:\n",
    "  - explained_fraction is the share of the observed infected-vs-baseline increase that lung-RNA\n",
    "    admixture alone is expected to contribute. Values well below 1 mean the increase is larger than\n",
    "    admixture can account for; values near or above 1 mean it cannot be distinguished from admixture.\n",
    "  - Inconsistent per-marker fractions mean a single admixture fraction describes the data poorly and\n",
    "    the model should be treated as qualitative only.\n",
    "  - This model cannot establish the SOURCE of the lung-derived signal.\n", sep = "")

cat("\nOutputs written to advanced_analyses/tables (tonsil_lung_admixture_*) and advanced_analyses/figures.\n")
