###############################################################################
# WGCNA candidate signature membership table
#
# Goal:
#   Cross-reference manuscript candidate genes with WGCNA module membership,
#   module-trait correlations, and intramodular kME.
#
# Why:
#   The generic top-hub tables are useful but not sufficient. We need to know
#   whether our biologically important candidates are strong members of the
#   infection/progression-associated modules.
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "openxlsx", "ggplot2", "pheatmap", "RColorBrewer")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}

library(tidyverse)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

candidate_file <- file.path(table_dir, "candidate_biomarker_table_full.csv")
if (!file.exists(candidate_file)) {
  stop("Run 01_candidate_biomarker_table.R first.")
}

candidates <- readr::read_csv(candidate_file, show_col_types = FALSE) %>%
  mutate(gene = as.character(gene)) %>%
  arrange(desc(priority_score), desc(n_significant_up)) %>%
  distinct(gene, .keep_all = TRUE)

load_wgcna_membership <- function(tissue_name) {
  rds_path <- file.path(table_dir, paste0("WGCNA_", tissue_name, "_object.rds"))
  gene_module_path <- file.path(table_dir, paste0("WGCNA_", tissue_name, "_gene_modules.csv"))
  corr_path <- file.path(table_dir, paste0("WGCNA_", tissue_name, "_module_trait_correlations.csv"))

  if (!file.exists(rds_path)) stop("Missing WGCNA object: ", rds_path)
  if (!file.exists(gene_module_path)) stop("Missing gene module table: ", gene_module_path)
  if (!file.exists(corr_path)) stop("Missing module-trait table: ", corr_path)

  obj <- readRDS(rds_path)
  datExpr <- obj$datExpr
  MEs <- obj$MEs

  gene_modules <- readr::read_csv(gene_module_path, show_col_types = FALSE) %>%
    mutate(gene = as.character(gene))

  candidate_genes_present <- intersect(candidates$gene, colnames(datExpr))

  kme_long <- map_dfr(candidate_genes_present, function(g) {
    module_num <- gene_modules %>%
      filter(gene == g) %>%
      pull(module_numeric) %>%
      .[1]
    module_col <- gene_modules %>%
      filter(gene == g) %>%
      pull(module_color) %>%
      .[1]
    me_name <- paste0("ME", module_num)
    if (!me_name %in% colnames(MEs)) return(NULL)

    tibble(
      tissue = tissue_name,
      gene = g,
      module_color = module_col,
      module_numeric = module_num,
      module_eigengene = me_name,
      kME = as.numeric(cor(datExpr[, g], MEs[, me_name], use = "p"))
    )
  })

  module_traits <- readr::read_csv(corr_path, show_col_types = FALSE) %>%
    filter(trait %in% c("infected", "dpi_numeric", "baseline", "dpi3", "dpi4", "dpi5")) %>%
    mutate(module_numeric = as.integer(str_remove(module, "^ME")))

  kme_long %>%
    left_join(
      module_traits %>%
        select(module_numeric, trait, correlation, pvalue, soft_power),
      by = "module_numeric"
    ) %>%
    mutate(tissue = tissue_name)
}

membership_long <- bind_rows(
  load_wgcna_membership("Lung"),
  load_wgcna_membership("Tonsil")
) %>%
  left_join(
    candidates %>%
      select(
        gene, category, manuscript_module, evidence_strength,
        n_significant_up, n_huvec_significant_up, n_in_vivo_significant_up,
        priority_score, conserved_call
      ),
    by = "gene"
  ) %>%
  mutate(
    candidate_is_strong_module_member = abs(kME) >= 0.8,
    module_trait_support = case_when(
      trait == "infected" & correlation >= 0.7 & pvalue < 0.01 ~ "strong_infection_positive",
      trait == "dpi_numeric" & correlation >= 0.7 & pvalue < 0.01 ~ "strong_progression_positive",
      trait == "infected" & correlation <= -0.7 & pvalue < 0.01 ~ "strong_infection_negative",
      trait == "dpi_numeric" & correlation <= -0.7 & pvalue < 0.01 ~ "strong_progression_negative",
      pvalue < 0.05 ~ "nominal_trait_association",
      TRUE ~ "weak_or_no_trait_association"
    )
  ) %>%
  arrange(tissue, module_numeric, desc(abs(kME)), gene, trait)

membership_summary <- membership_long %>%
  filter(trait %in% c("infected", "dpi_numeric")) %>%
  group_by(
    tissue, gene, category, manuscript_module, evidence_strength,
    module_color, module_numeric, module_eigengene, kME,
    n_significant_up, n_huvec_significant_up, n_in_vivo_significant_up,
    priority_score, conserved_call
  ) %>%
  summarise(
    infected_module_cor = correlation[trait == "infected"][1],
    infected_module_p = pvalue[trait == "infected"][1],
    dpi_module_cor = correlation[trait == "dpi_numeric"][1],
    dpi_module_p = pvalue[trait == "dpi_numeric"][1],
    strong_module_member = abs(kME) >= 0.8,
    strong_infection_module = !is.na(infected_module_cor) &
      infected_module_cor >= 0.7 & infected_module_p < 0.01,
    strong_progression_module = !is.na(dpi_module_cor) &
      dpi_module_cor >= 0.7 & dpi_module_p < 0.01,
    .groups = "drop"
  ) %>%
  mutate(
    wgcna_support_class = case_when(
      strong_module_member & strong_infection_module & strong_progression_module ~
        "strong infection/progression module member",
      strong_module_member & strong_infection_module ~
        "strong infection module member",
      strong_module_member & strong_progression_module ~
        "strong progression module member",
      strong_module_member ~
        "strong module member but weak trait support",
      TRUE ~ "present but not strong module member"
    )
  ) %>%
  arrange(tissue, desc(strong_infection_module), desc(abs(kME)), gene)

membership_summary <- membership_summary %>%
  distinct(tissue, gene, module_numeric, module_eigengene, kME, .keep_all = TRUE)

readr::write_csv(
  membership_long,
  file.path(table_dir, "WGCNA_candidate_signature_membership_long.csv")
)

readr::write_csv(
  membership_summary,
  file.path(table_dir, "WGCNA_candidate_signature_membership_summary.csv")
)

wb <- createWorkbook()
addWorksheet(wb, "summary")
writeData(wb, "summary", membership_summary)
addWorksheet(wb, "long")
writeData(wb, "long", membership_long)
freezePane(wb, "summary", firstRow = TRUE)
freezePane(wb, "long", firstRow = TRUE)
saveWorkbook(
  wb,
  file.path(table_dir, "WGCNA_candidate_signature_membership.xlsx"),
  overwrite = TRUE
)

###############################################################################
# Figures
###############################################################################

plot_df <- membership_summary %>%
  mutate(
    gene = fct_reorder(gene, abs(kME)),
    tissue = factor(tissue, levels = c("Lung", "Tonsil")),
    manuscript_module = factor(
      manuscript_module,
      levels = c(
        "Conserved antiviral/IFN core",
        "In vivo complement/coagulation disease module",
        "Mainly in vivo immune/progression module",
        "HUVEC-enriched endothelial/early-response module",
        "Supporting candidate"
      )
    )
  )

png(
  file.path(figure_dir, "WGCNA_candidate_signature_membership_dotplot.png"),
  width = 3400, height = 2600, res = 300
)
print(
  ggplot(plot_df, aes(x = tissue, y = gene)) +
    geom_point(
      aes(size = abs(kME), fill = infected_module_cor, shape = strong_module_member),
      color = "grey20",
      stroke = 0.25,
      alpha = 0.95
    ) +
    facet_grid(manuscript_module ~ ., scales = "free_y", space = "free_y") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24)) +
    scale_size_continuous(range = c(1.5, 7), limits = c(0, 1), name = "|kME|") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Module-infection\ncorrelation"
    ) +
    labs(
      title = "WGCNA support for candidate Nipah host-response genes",
      subtitle = "Point size shows module membership strength; color shows module correlation with infection",
      x = NULL,
      y = NULL,
      shape = "Strong module\nmember"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text.y = element_text(angle = 0, face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
)
dev.off()

heatmap_df <- membership_summary %>%
  mutate(tissue_gene = paste(tissue, gene, sep = ": ")) %>%
  select(tissue_gene, infected_module_cor, dpi_module_cor, kME) %>%
  column_to_rownames("tissue_gene") %>%
  as.matrix()

png(
  file.path(figure_dir, "WGCNA_candidate_signature_membership_heatmap.png"),
  width = 2200, height = 2800, res = 300
)
pheatmap(
  heatmap_df,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101),
  breaks = seq(-1, 1, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize = 7,
  main = "Candidate gene WGCNA membership and module-trait support"
)
dev.off()

message("WGCNA candidate membership analysis complete.")
message("Summary table: ", file.path(table_dir, "WGCNA_candidate_signature_membership_summary.csv"))
