###############################################################################
# GSE310471 ordered-time (trend) analysis across days post-infection (DPI)
#
# Rationale:
#   The DESeq2 LRT (script 18) treats DPI as a categorical factor and therefore
#   tests whether expression differs among time points. It does not show a
#   monotonic temporal trend. This script adds formal ordered-time analyses:
#
#   (A) Module eigengenes (lung tan/ME12; tonsil blue/ME2; tonsil ME2 recomputed
#       without lung genes) versus numeric DPI:
#         - Spearman correlation, all samples (baseline = 0, 3, 4, 5 DPI)
#         - Spearman correlation, infected samples only (3, 4, 5 DPI), which
#           separates an infected-versus-baseline step from a graded DPI trend
#         - Jonckheere-Terpstra test for an ordered alternative
#           (baseline < 3 < 4 < 5 DPI), permutation P value
#         - Linear slope per DPI (all samples, infected only) with 95% CI
#         - Wilcoxon test, infected versus baseline
#   (B) Gene-level numeric-DPI Wald test (DESeq2, design ~ dpi_numeric) per
#       tissue, in all samples and in infected samples only, with the number
#       of genes passing padj < 0.05 and their direction.
#
# DPI is a time-since-inoculation variable. It is not disease severity: no
# pathology, clinical score or viral-load covariates are available in the
# public metadata analysed here.
#
# Inputs:
#   GSE310471/results/GSE310471_analysis_objects.rds
#   GSE310471/results/tables/GSE310471_sample_metadata.csv
#   advanced_analyses/tables/WGCNA_Lung_object.rds, WGCNA_Tonsil_object.rds
#   advanced_analyses/tables/tonsil_ME2_per_sample_eigengenes.csv (script 22)
#
# Outputs (advanced_analyses/tables, advanced_analyses/figures):
#   GSE310471_ordered_time_module_trend_tests.csv
#   GSE310471_ordered_time_gene_level_summary.csv
#   GSE310471_ordered_time_gene_level_full.csv
#   GSE310471_ordered_time_module_trend.png
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(20261003)

packages <- c("tidyverse", "DESeq2", "ggplot2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "DESeq2") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({ library(tidyverse); library(DESeq2); library(ggplot2) })

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
gse_dir <- file.path(project_dir, "GSE310471")
table_dir <- file.path(project_dir, "advanced_analyses", "tables")
figure_dir <- file.path(project_dir, "advanced_analyses", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- readr::read_csv(file.path(gse_dir, "results", "tables", "GSE310471_sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(dpi_numeric = ifelse(dpi == "baseline", 0, as.numeric(str_extract(dpi, "[0-9]+"))),
         infected = as.integer(dpi != "baseline"))
meta_sel <- metadata %>% select(sample_id, dpi, dpi_numeric, infected)  # drop metadata 'tissue' to avoid join suffixes

###############################################################################
# (A) Module eigengene trend tests
###############################################################################

# Jonckheere-Terpstra statistic (ordered groups), permutation P value
jt_stat <- function(x, g) {
  lv <- sort(unique(g)); s <- 0
  for (i in seq_len(length(lv) - 1)) for (j in (i + 1):length(lv)) {
    a <- x[g == lv[i]]; b <- x[g == lv[j]]
    s <- s + sum(outer(a, b, function(u, v) (v > u) + 0.5 * (v == u)))
  }
  s
}
jt_test <- function(x, g, nperm = 20000) {
  obs <- jt_stat(x, g)
  perm <- replicate(nperm, jt_stat(x, sample(g)))
  p_up <- (sum(perm >= obs) + 1) / (nperm + 1)
  p_dn <- (sum(perm <= obs) + 1) / (nperm + 1)
  tibble(JT_stat = obs, JT_p_increasing = p_up, JT_p_two_sided = min(1, 2 * min(p_up, p_dn)))
}
slope_ci <- function(y, x) {
  if (length(unique(x)) < 2) return(tibble(slope = NA, lo = NA, hi = NA))
  f <- lm(y ~ x); ci <- suppressWarnings(confint(f))["x", ]
  tibble(slope = unname(coef(f)["x"]), lo = ci[1], hi = ci[2])
}
sp <- function(x, y) { r <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)); c(rho = unname(r$estimate), p = r$p.value) }

load_me <- function(rds, me_col, tissue, label) {
  obj <- readRDS(rds); MEs <- obj$MEs
  tibble(sample_id = rownames(MEs), eigengene = MEs[[me_col]], tissue = tissue, module = label) %>%
    left_join(meta_sel, by = "sample_id")
}
lung_me <- load_me(file.path(table_dir, "WGCNA_Lung_object.rds"), "ME12", "Lung", "Lung tan/ME12 (antiviral module)")
tonsil_me <- load_me(file.path(table_dir, "WGCNA_Tonsil_object.rds"), "ME2", "Tonsil", "Tonsil blue/ME2 (complement/coagulation module)")

me_all <- bind_rows(lung_me, tonsil_me)
ex_path <- file.path(table_dir, "tonsil_ME2_per_sample_eigengenes.csv")
if (file.exists(ex_path)) {
  ex <- readr::read_csv(ex_path, show_col_types = FALSE)
  me_all <- bind_rows(me_all,
    ex %>% transmute(sample_id, eigengene = ME2_lung_excluded, tissue = "Tonsil",
                     module = "Tonsil blue/ME2, lung genes excluded") %>% left_join(meta_sel, by = "sample_id"))
}
stopifnot(!anyNA(me_all$dpi_numeric))

trend_tbl <- me_all %>% group_by(module, tissue) %>% group_modify(~{
  d <- .x; di <- d %>% filter(infected == 1)
  sa <- sp(d$dpi_numeric, d$eigengene); si <- sp(di$dpi_numeric, di$eigengene)
  jt <- jt_test(d$eigengene, d$dpi_numeric)
  jti <- jt_test(di$eigengene, di$dpi_numeric)
  sl <- slope_ci(d$eigengene, d$dpi_numeric); sli <- slope_ci(di$eigengene, di$dpi_numeric)
  wt <- wilcox.test(eigengene ~ infected, data = d, exact = FALSE)
  tibble(
    n_all = nrow(d), n_infected = nrow(di),
    spearman_rho_all = sa["rho"], spearman_p_all = sa["p"],
    spearman_rho_infected_only = si["rho"], spearman_p_infected_only = si["p"],
    JT_p_increasing_all = jt$JT_p_increasing, JT_p_two_sided_all = jt$JT_p_two_sided,
    JT_p_increasing_infected_only = jti$JT_p_increasing, JT_p_two_sided_infected_only = jti$JT_p_two_sided,
    slope_per_DPI_all = sl$slope, slope_all_CI_low = sl$lo, slope_all_CI_high = sl$hi,
    slope_per_DPI_infected_only = sli$slope, slope_infected_CI_low = sli$lo, slope_infected_CI_high = sli$hi,
    median_baseline = median(d$eigengene[d$infected == 0]), median_infected = median(di$eigengene),
    wilcoxon_p_infected_vs_baseline = wt$p.value
  )
}) %>% ungroup()
readr::write_csv(trend_tbl, file.path(table_dir, "GSE310471_ordered_time_module_trend_tests.csv"))

###############################################################################
# (B) Gene-level numeric-DPI Wald test (DESeq2)
###############################################################################

analysis_objects <- readRDS(file.path(gse_dir, "results", "GSE310471_analysis_objects.rds"))
run_numeric <- function(tissue_name, infected_only) {
  dds0 <- analysis_objects$tissue_results[[tissue_name]]$dds
  cd <- as.data.frame(colData(dds0)) %>% tibble::rownames_to_column("sample_id") %>%
    select(sample_id) %>% left_join(meta_sel, by = "sample_id")
  keep <- if (infected_only) cd$infected == 1 else rep(TRUE, nrow(cd))
  cnt <- counts(dds0)[, keep, drop = FALSE]
  coldat <- data.frame(dpi_numeric = cd$dpi_numeric[keep], row.names = cd$sample_id[keep])
  d <- DESeqDataSetFromMatrix(countData = round(cnt), colData = coldat, design = ~ dpi_numeric)
  d <- d[rowSums(counts(d)) >= 10, ]
  d <- DESeq(d, quiet = TRUE)
  res <- results(d, name = "dpi_numeric") %>% as.data.frame() %>% tibble::rownames_to_column("Gene")
  res$tissue <- tissue_name
  res$samples <- ifelse(infected_only, "infected only (3-5 DPI)", "all (baseline + 3-5 DPI)")
  res$n_samples <- ncol(d)
  res
}
gene_full <- bind_rows(lapply(c("Lung", "Tonsil"), function(t)
  bind_rows(run_numeric(t, FALSE), run_numeric(t, TRUE))))
readr::write_csv(gene_full, file.path(table_dir, "GSE310471_ordered_time_gene_level_full.csv"))

gene_sum <- gene_full %>% group_by(tissue, samples, n_samples) %>%
  summarise(n_genes_tested = sum(!is.na(padj)),
            n_sig_padj_0.05 = sum(padj < 0.05, na.rm = TRUE),
            n_sig_increasing = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
            n_sig_decreasing = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE),
            pct_sig = round(100 * n_sig_padj_0.05 / n_genes_tested, 2), .groups = "drop")
readr::write_csv(gene_sum, file.path(table_dir, "GSE310471_ordered_time_gene_level_summary.csv"))

###############################################################################
# Figure: eigengene versus numeric DPI with infected-only linear fit
###############################################################################

p <- ggplot(me_all, aes(x = dpi_numeric, y = eigengene)) +
  geom_jitter(aes(colour = factor(infected, labels = c("Baseline", "Infected"))), width = 0.12, size = 2, alpha = 0.85) +
  stat_summary(aes(group = dpi_numeric), fun = median, geom = "crossbar", width = 0.35, colour = "grey30", linewidth = 0.3) +
  geom_smooth(data = ~ filter(.x, infected == 1), method = "lm", se = TRUE, colour = "black", linewidth = 0.6, alpha = 0.15) +
  facet_wrap(~module, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = c(0, 3, 4, 5), labels = c("baseline", "3", "4", "5")) +
  scale_colour_manual(values = c(Baseline = "#6c757d", Infected = "#1f4e79"), name = NULL) +
  labs(title = "Module eigengenes versus infection time",
       subtitle = "Line: linear fit across infected samples only (3-5 DPI); bars: group medians",
       x = "Days post-infection", y = "Module eigengene") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "GSE310471_ordered_time_module_trend.png"), p, width = 7, height = 3 + 2.4 * length(unique(me_all$module)), dpi = 300)

cat("\n=== Module trend tests ===\n"); print(as.data.frame(trend_tbl %>% mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)
cat("\n=== Gene-level numeric-DPI summary ===\n"); print(as.data.frame(gene_sum), row.names = FALSE)
