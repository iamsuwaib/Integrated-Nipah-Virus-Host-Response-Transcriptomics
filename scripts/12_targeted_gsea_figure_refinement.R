###############################################################################
# 08. Targeted GSEA and WGCNA enrichment figure refinement
###############################################################################

options(stringsAsFactors = FALSE)

cran_pkgs <- c(
  "tidyverse", "ggplot2", "readr", "dplyr",
  "stringr", "forcats", "scales", "grid"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(cran_pkgs, install_if_missing))

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
})

# Avoid namespace conflicts
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange

# ----------------------------- Paths -----------------------------------------

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")

gsea_table_dir <- file.path(advanced_dir, "tables", "gsea_module_enrichment")
out_figure_dir <- file.path(advanced_dir, "figures", "gsea_module_enrichment", "targeted_refined")
out_table_dir <- file.path(advanced_dir, "tables", "gsea_module_enrichment", "targeted_refined")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

gsea_file <- file.path(gsea_table_dir, "GO_BP_GSEA_all_contrasts.csv")
module_file <- file.path(gsea_table_dir, "GO_BP_WGCNA_module_ORA_all.csv")

if (!file.exists(gsea_file)) stop("Missing GSEA table: ", gsea_file)
if (!file.exists(module_file)) stop("Missing WGCNA module ORA table: ", module_file)

gsea_all <- readr::read_csv(gsea_file, show_col_types = FALSE)
module_all <- readr::read_csv(module_file, show_col_types = FALSE)

# ----------------------------- Settings --------------------------------------

bio_theme_order <- c(
  "Antiviral/IFN",
  "Complement",
  "Coagulation/hemostasis",
  "Vascular/leukocyte",
  "Metabolic remodeling"
)

bio_theme_short_labels <- c(
  "Antiviral/IFN" = "Antiviral\nIFN",
  "Complement" = "Complement",
  "Coagulation/hemostasis" = "Coagulation\nhemostasis",
  "Vascular/leukocyte" = "Vascular\nleukocyte",
  "Metabolic remodeling" = "Metabolic\nremodeling"
)

# ------------------------- Target GO terms -----------------------------------

target_terms <- tibble::tribble(
  ~bio_theme, ~Description, ~display_term,
  "Antiviral/IFN", "defense response to virus", "Defense response to virus",
  "Antiviral/IFN", "response to virus", "Response to virus",
  "Antiviral/IFN", "response to type I interferon", "Response to type I interferon",
  "Antiviral/IFN", "type I interferon-mediated signaling pathway", "Type I IFN signaling",
  "Antiviral/IFN", "interferon-mediated signaling pathway", "Interferon signaling",
  "Antiviral/IFN", "antiviral innate immune response", "Antiviral innate immune response",
  
  "Complement", "complement activation", "Complement activation",
  "Complement", "complement activation, classical pathway", "Complement, classical pathway",
  "Complement", "complement activation, alternative pathway", "Complement, alternative pathway",
  
  "Coagulation/hemostasis", "regulation of blood coagulation", "Regulation of blood coagulation",
  "Coagulation/hemostasis", "blood coagulation", "Blood coagulation",
  "Coagulation/hemostasis", "hemostasis", "Hemostasis",
  "Coagulation/hemostasis", "fibrinolysis", "Fibrinolysis",
  
  "Vascular/leukocyte", "vascular process in circulatory system", "Vascular process",
  "Vascular/leukocyte", "leukocyte migration", "Leukocyte migration",
  "Vascular/leukocyte", "endothelial cell migration", "Endothelial cell migration",
  
  "Metabolic remodeling", "oxidative phosphorylation", "Oxidative phosphorylation",
  "Metabolic remodeling", "mitochondrial ATP synthesis coupled electron transport", "Mitochondrial ATP synthesis",
  "Metabolic remodeling", "aerobic respiration", "Aerobic respiration",
  "Metabolic remodeling", "ribosome biogenesis", "Ribosome biogenesis"
) %>%
  dplyr::mutate(
    bio_theme = factor(bio_theme, levels = bio_theme_order),
    Description_key = stringr::str_squish(stringr::str_to_lower(Description))
  )

contrast_labels <- tibble::tribble(
  ~dataset, ~contrast, ~contrast_label, ~contrast_group, ~contrast_order,
  "GSE32902", "NiV_vs_Mock", "HUVEC\nNiV/mock", "HUVEC", 1,
  "GSE33133", "NiV_vs_Mock", "HUVEC-ext\nNiV/mock", "HUVEC", 2,
  "GSE33133", "NiVdC_vs_Mock", "HUVEC-ext\nNiV-dC/mock", "HUVEC", 3,
  "GSE33133", "NiVdC_vs_NiV", "HUVEC-ext\nNiV-dC/NiV", "HUVEC", 4,
  "GSE310471", "Lung_3DPI_vs_baseline", "Lung\n3 DPI", "Lung", 5,
  "GSE310471", "Lung_4DPI_vs_baseline", "Lung\n4 DPI", "Lung", 6,
  "GSE310471", "Lung_5DPI_vs_baseline", "Lung\n5 DPI", "Lung", 7,
  "GSE310471", "Tonsil_3DPI_vs_baseline", "Tonsil\n3 DPI", "Tonsil", 8,
  "GSE310471", "Tonsil_4DPI_vs_baseline", "Tonsil\n4 DPI", "Tonsil", 9,
  "GSE310471", "Tonsil_5DPI_vs_baseline", "Tonsil\n5 DPI", "Tonsil", 10
)

# ------------------------- Targeted GSEA -------------------------------------

gsea_targeted <- gsea_all %>%
  dplyr::mutate(
    Description_key = stringr::str_squish(stringr::str_to_lower(Description))
  ) %>%
  dplyr::inner_join(target_terms, by = "Description_key", suffix = c("", "_target")) %>%
  dplyr::inner_join(contrast_labels, by = c("dataset", "contrast")) %>%
  dplyr::mutate(
    padj_plot = pmax(p.adjust, 1e-300),
    neg_log10_fdr = -log10(padj_plot),
    fdr_label = dplyr::case_when(
      p.adjust <= 0.001 ~ "FDR <= 0.001",
      p.adjust <= 0.01 ~ "FDR <= 0.01",
      p.adjust <= 0.05 ~ "FDR <= 0.05",
      p.adjust <= 0.25 ~ "FDR <= 0.25",
      TRUE ~ "FDR > 0.25"
    ),
    display_term = factor(display_term, levels = rev(unique(target_terms$display_term))),
    bio_theme = factor(bio_theme, levels = bio_theme_order),
    contrast_label = factor(contrast_label, levels = contrast_labels$contrast_label),
    contrast_group = factor(contrast_group, levels = c("HUVEC", "Lung", "Tonsil"))
  )

readr::write_csv(
  gsea_targeted %>%
    dplyr::arrange(bio_theme, display_term, contrast_order) %>%
    dplyr::select(
      dataset, system, contrast, contrast_label, contrast_group,
      bio_theme, display_term, ID, Description, NES, pvalue,
      p.adjust, qvalue, setSize, direction
    ),
  file.path(out_table_dir, "targeted_GO_BP_GSEA_terms_for_figure.csv")
)

gsea_plot_data <- gsea_targeted %>%
  dplyr::filter(!is.na(p.adjust), p.adjust <= 0.25) %>%
  dplyr::mutate(
    neg_log10_fdr_capped = pmin(neg_log10_fdr, 25),
    NES_capped = pmax(pmin(NES, 3.2), -3.2)
  )

if (nrow(gsea_plot_data) > 0) {
  
  gsea_plot <- ggplot(
    gsea_plot_data,
    aes(
      x = contrast_label,
      y = display_term,
      size = neg_log10_fdr_capped,
      color = NES_capped
    )
  ) +
    geom_point(alpha = 0.92) +
    facet_grid(
      bio_theme ~ contrast_group,
      scales = "free",
      space = "free",
      switch = "y",
      labeller = labeller(bio_theme = bio_theme_short_labels)
    ) +
    scale_color_gradient2(
      low = "#2B6CB0",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-3.2, 3.2),
      oob = scales::squish,
      name = "NES"
    ) +
    scale_size_continuous(
      range = c(2.2, 8.5),
      breaks = c(2, 5, 10, 20),
      name = expression(-log[10]("FDR"))
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9.5) +
    theme(
      text = element_text(color = "black"),
      strip.background = element_rect(fill = "#F1F3F5", color = "#BFC5CC", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = 9.5),
      strip.text.y.left = element_text(
        angle = 0,
        hjust = 0.5,
        lineheight = 0.95,
        margin = margin(2, 3, 2, 3)
      ),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.2),
      axis.text.y = element_text(size = 8.7),
      panel.grid.major = element_line(color = "#E6E8EB", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.45, "lines"),
      panel.spacing.y = unit(0.32, "lines"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8.8),
      legend.text = element_text(size = 8),
      legend.key.height = unit(5, "mm"),
      legend.key.width = unit(5, "mm"),
      plot.margin = margin(8, 10, 8, 8)
    ) +
    guides(
      color = guide_colorbar(order = 1, barheight = unit(45, "mm")),
      size = guide_legend(order = 2, override.aes = list(color = "#4A5568"))
    )
  
  ggsave(
    filename = file.path(out_figure_dir, "Figure_targeted_GO_BP_GSEA_dotplot.png"),
    plot = gsea_plot,
    width = 10.8,
    height = 7.4,
    dpi = 600,
    bg = "white"
  )
  
} else {
  warning("No targeted GSEA terms with p.adjust <= 0.25. GSEA dotplot was not generated.")
}

# ------------------------- Targeted WGCNA ORA --------------------------------

module_target_terms <- tibble::tribble(
  ~module_focus, ~Description, ~display_term,
  "Lung antiviral module", "defense response to virus", "Defense response to virus",
  "Lung antiviral module", "response to virus", "Response to virus",
  "Lung antiviral module", "response to type I interferon", "Response to type I interferon",
  "Lung antiviral module", "type I interferon-mediated signaling pathway", "Type I IFN signaling",
  "Lung antiviral module", "interferon-mediated signaling pathway", "Interferon signaling",
  
  "Tonsil progression module", "innate immune response", "Innate immune response",
  "Tonsil progression module", "complement activation", "Complement activation",
  "Tonsil progression module", "complement activation, classical pathway", "Complement, classical pathway",
  "Tonsil progression module", "complement activation, alternative pathway", "Complement, alternative pathway",
  "Tonsil progression module", "blood coagulation", "Blood coagulation",
  "Tonsil progression module", "hemostasis", "Hemostasis",
  "Tonsil progression module", "regulation of blood coagulation", "Regulation of blood coagulation",
  "Tonsil progression module", "negative regulation of hemostasis", "Negative regulation of hemostasis",
  "Tonsil progression module", "fibrinolysis", "Fibrinolysis"
) %>%
  dplyr::mutate(
    Description_key = stringr::str_squish(stringr::str_to_lower(Description))
  )

module_targeted <- module_all %>%
  dplyr::mutate(
    Description_key = stringr::str_squish(stringr::str_to_lower(Description))
  ) %>%
  dplyr::inner_join(module_target_terms, by = "Description_key", suffix = c("", "_target")) %>%
  dplyr::mutate(
    module_panel = dplyr::case_when(
      tissue == "Lung" ~ "Lung tan/ME12: antiviral IFN module",
      tissue == "Tonsil" ~ "Tonsil blue/ME2: progression module",
      TRUE ~ module
    ),
    padj_plot = pmax(p.adjust, 1e-300),
    neg_log10_fdr = -log10(padj_plot),
    display_term = factor(display_term, levels = rev(unique(module_target_terms$display_term))),
    module_panel = factor(
      module_panel,
      levels = c(
        "Lung tan/ME12: antiviral IFN module",
        "Tonsil blue/ME2: progression module"
      )
    )
  )

readr::write_csv(
  module_targeted %>%
    dplyr::arrange(tissue, p.adjust) %>%
    dplyr::select(
      tissue, module, module_panel, display_term,
      ID, Description, GeneRatio, BgRatio, Count,
      pvalue, p.adjust, qvalue, geneID,
      module_genes_input, module_genes_mapped,
      universe_genes_input, universe_genes_mapped
    ),
  file.path(out_table_dir, "targeted_WGCNA_module_ORA_terms_for_figure.csv")
)

module_plot_data <- module_targeted %>%
  dplyr::filter(!is.na(p.adjust), p.adjust <= 0.25) %>%
  dplyr::mutate(
    neg_log10_fdr_capped = pmin(neg_log10_fdr, 25),
    display_term = forcats::fct_reorder(display_term, neg_log10_fdr, .desc = FALSE)
  )

if (nrow(module_plot_data) > 0) {
  
  module_plot <- ggplot(
    module_plot_data,
    aes(x = neg_log10_fdr_capped, y = display_term)
  ) +
    geom_col(aes(fill = tissue), width = 0.68, alpha = 0.92) +
    geom_text(
      aes(label = Count),
      hjust = -0.15,
      size = 3.5,
      color = "black"
    ) +
    facet_wrap(~ module_panel, scales = "free_y", ncol = 1) +
    scale_fill_manual(
      values = c("Lung" = "#2B6CB0", "Tonsil" = "#B2182B"),
      guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = expression(-log[10]("FDR")), y = NULL) +
    theme_bw(base_size = 9.5) +
    theme(
      text = element_text(color = "black"),
      strip.background = element_rect(fill = "#F1F3F5", color = "#BFC5CC", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = 9.5),
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 8.7),
      axis.title.x = element_text(face = "bold", size = 9.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0.55, "lines"),
      plot.margin = margin(8, 18, 8, 8)
    )
  
  ggsave(
    filename = file.path(out_figure_dir, "Figure_targeted_WGCNA_module_ORA_barplot.png"),
    plot = module_plot,
    width = 7.2,
    height = 6.4,
    dpi = 600,
    bg = "white"
  )
  
} else {
  warning("No targeted WGCNA module ORA terms with p.adjust <= 0.25. Module barplot was not generated.")
}

message("Targeted GSEA and WGCNA enrichment figure refinement completed.")
message("Figures written to: ", out_figure_dir)
message("Tables written to: ", out_table_dir)