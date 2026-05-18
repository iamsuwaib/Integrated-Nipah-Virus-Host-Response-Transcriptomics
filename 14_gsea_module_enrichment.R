###############################################################################
# 14. GSEA and WGCNA module enrichment for Nipah host-response manuscript
#
# Purpose:
#   1) Run ranked GO Biological Process GSEA for each analyzed contrast.
#   2) Run GO Biological Process over-representation analysis for prioritized
#      WGCNA modules: lung tan/ME12 and tonsil blue/ME2.
#   3) Export clean tables and journal-friendly PNG summary figures.
#
# Notes:
#   - This script does not pool raw expression across datasets.
#   - Each contrast is ranked from its own differential-expression result table.
#   - For microarray datasets, moderated t statistics are used as the primary
#     ranking statistic.
#   - For RNA-seq contrasts, DESeq2 Wald statistics are used as the primary
#     ranking statistic.
#   - WGCNA module enrichment is treated as module-level biological annotation,
#     not as mechanistic proof.
###############################################################################

options(stringsAsFactors = FALSE)

# ------------------------------ 1. Packages ----------------------------------

cran_pkgs <- c(
  "tidyverse", "ggplot2", "readr", "dplyr", "stringr", "forcats",
  "patchwork", "scales"
)
bioc_pkgs <- c(
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% bioc_pkgs) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

invisible(lapply(c(cran_pkgs, bioc_pkgs), install_if_missing))

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(patchwork)
  library(scales)
})

# ------------------------------ 2. Paths -------------------------------------

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")

out_table_dir <- file.path(table_dir, "gsea_module_enrichment")
out_figure_dir <- file.path(figure_dir, "gsea_module_enrichment")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- 3. Settings -----------------------------------

set.seed(20260517)

GSEA_MIN_SIZE <- 10
GSEA_MAX_SIZE <- 500
GSEA_EPS <- 1e-10
GSEA_P_CUTOFF <- 1
GSEA_PADJ_REPORT <- 0.25

ORA_MIN_SIZE <- 10
ORA_MAX_SIZE <- 500
ORA_P_CUTOFF <- 1
ORA_PADJ_REPORT <- 0.25

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey25", size = base_size),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      strip.text = element_text(face = "bold")
    )
}

clean_term <- function(x, max_chars = 58) {
  x %>%
    str_replace_all("_", " ") %>%
    str_to_sentence() %>%
    str_wrap(width = max_chars)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

collapse_ranked_genes <- function(df, entrez_col, stat_col, fc_col = NULL) {
  entrez <- df[[entrez_col]]
  stat <- df[[stat_col]]

  ranked <- tibble(
    ENTREZID = as.character(entrez),
    rank_stat = as.numeric(stat),
    fallback_fc = if (!is.null(fc_col) && fc_col %in% names(df)) {
      as.numeric(df[[fc_col]])
    } else {
      NA_real_
    }
  ) %>%
    filter(!is.na(ENTREZID), ENTREZID != "", ENTREZID != "NA") %>%
    mutate(
      ENTREZID = str_split(ENTREZID, pattern = "[,;/ ]+", simplify = TRUE)[, 1],
      rank_stat = if_else(is.na(rank_stat), fallback_fc, rank_stat)
    ) %>%
    filter(!is.na(rank_stat), is.finite(rank_stat)) %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(rank_stat))

  gene_list <- ranked$rank_stat
  names(gene_list) <- ranked$ENTREZID
  sort(gene_list, decreasing = TRUE)
}

map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) return(tibble(SYMBOL = character(), ENTREZID = character()))

  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    as_tibble() %>%
    filter(!is.na(ENTREZID), ENTREZID != "") %>%
    distinct(SYMBOL, ENTREZID)
}

run_gsea_go <- function(gene_list, contrast_label, dataset_label, system_label) {
  message("Running GO-GSEA: ", contrast_label, " (", length(gene_list), " ranked genes)")

  if (length(gene_list) < 100) {
    warning("Skipping ", contrast_label, ": fewer than 100 ranked Entrez genes.")
    return(tibble())
  }

  ego <- tryCatch(
    {
      clusterProfiler::gseGO(
        geneList = gene_list,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        keyType = "ENTREZID",
        minGSSize = GSEA_MIN_SIZE,
        maxGSSize = GSEA_MAX_SIZE,
        eps = GSEA_EPS,
        pvalueCutoff = GSEA_P_CUTOFF,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
    },
    error = function(e) {
      warning("GSEA failed for ", contrast_label, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(tibble())

  as.data.frame(ego) %>%
    as_tibble() %>%
    mutate(
      dataset = dataset_label,
      system = system_label,
      contrast = contrast_label,
      direction = if_else(NES >= 0, "Enriched among upregulated genes", "Enriched among downregulated genes"),
      term_label = clean_term(Description)
    ) %>%
    select(dataset, system, contrast, ID, Description, term_label, setSize,
           enrichmentScore, NES, pvalue, p.adjust, qvalue, rank,
           leading_edge, core_enrichment, direction)
}

run_module_ora_go <- function(module_genes, universe_genes, module_label, tissue_label) {
  message("Running module ORA: ", module_label)

  module_map <- map_symbols_to_entrez(module_genes)
  universe_map <- map_symbols_to_entrez(universe_genes)

  module_entrez <- unique(module_map$ENTREZID)
  universe_entrez <- unique(universe_map$ENTREZID)

  if (length(module_entrez) < 5) {
    warning("Skipping module ORA for ", module_label, ": fewer than 5 mapped genes.")
    return(tibble())
  }

  ego <- tryCatch(
    {
      clusterProfiler::enrichGO(
        gene = module_entrez,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = ORA_P_CUTOFF,
        qvalueCutoff = ORA_P_CUTOFF,
        minGSSize = ORA_MIN_SIZE,
        maxGSSize = ORA_MAX_SIZE,
        readable = TRUE
      )
    },
    error = function(e) {
      warning("Module ORA failed for ", module_label, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(tibble())

  as.data.frame(ego) %>%
    as_tibble() %>%
    mutate(
      tissue = tissue_label,
      module = module_label,
      term_label = clean_term(Description),
      module_genes_input = length(unique(module_genes)),
      module_genes_mapped = length(module_entrez),
      universe_genes_input = length(unique(universe_genes)),
      universe_genes_mapped = length(universe_entrez)
    ) %>%
    select(tissue, module, ID, Description, term_label, GeneRatio, BgRatio,
           pvalue, p.adjust, qvalue, geneID, Count,
           module_genes_input, module_genes_mapped,
           universe_genes_input, universe_genes_mapped)
}

# ------------------------- 4. Contrast-level GSEA ----------------------------

contrast_specs <- tribble(
  ~dataset, ~system, ~contrast, ~path, ~entrez_col, ~stat_col, ~fc_col,
  "GSE32902", "HUVEC endothelial layer", "NiV_vs_Mock",
  file.path(project_dir, "GSE32902/results/tables/GSE32902_limma_all_probes_annotated.csv"),
  "ENTREZID", "t", "logFC",

  "GSE33133", "HUVEC endothelial layer", "NiV_vs_Mock",
  file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiV_vs_Mock_annotated.csv"),
  "ENTREZID", "t", "logFC",

  "GSE33133", "HUVEC endothelial layer", "NiVdC_vs_Mock",
  file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_Mock_annotated.csv"),
  "ENTREZID", "t", "logFC",

  "GSE33133", "HUVEC endothelial layer", "NiVdC_vs_NiV",
  file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_NiV_annotated.csv"),
  "ENTREZID", "t", "logFC",

  "GSE310471", "AGM lung in vivo layer", "Lung_3DPI_vs_baseline",
  file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_3DPI_vs_baseline.csv"),
  "human_ENTREZID", "stat", "log2FoldChange",

  "GSE310471", "AGM lung in vivo layer", "Lung_4DPI_vs_baseline",
  file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_4DPI_vs_baseline.csv"),
  "human_ENTREZID", "stat", "log2FoldChange",

  "GSE310471", "AGM lung in vivo layer", "Lung_5DPI_vs_baseline",
  file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_5DPI_vs_baseline.csv"),
  "human_ENTREZID", "stat", "log2FoldChange",

  "GSE310471", "AGM tonsil in vivo layer", "Tonsil_3DPI_vs_baseline",
  file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_3DPI_vs_baseline.csv"),
  "human_ENTREZID", "stat", "log2FoldChange",

  "GSE310471", "AGM tonsil in vivo layer", "Tonsil_4DPI_vs_baseline",
  file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_4DPI_vs_baseline.csv"),
  "human_ENTREZID", "stat", "log2FoldChange",

  "GSE310471", "AGM tonsil in vivo layer", "Tonsil_5DPI_vs_baseline",
  file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_5DPI_vs_baseline.csv"),
  "human_ENTREZID", "stat", "log2FoldChange"
)

rank_summary <- list()
gsea_results <- pmap_dfr(
  contrast_specs,
  function(dataset, system, contrast, path, entrez_col, stat_col, fc_col) {
    df <- safe_read_csv(path)
    gene_list <- collapse_ranked_genes(df, entrez_col, stat_col, fc_col)

    rank_summary[[paste(dataset, contrast, sep = "__")]] <<- tibble(
      dataset = dataset,
      system = system,
      contrast = contrast,
      input_rows = nrow(df),
      ranked_entrez_genes = length(gene_list),
      positive_ranked_genes = sum(gene_list > 0),
      negative_ranked_genes = sum(gene_list < 0),
      max_rank_stat = max(gene_list, na.rm = TRUE),
      min_rank_stat = min(gene_list, na.rm = TRUE)
    )

    readr::write_csv(
      tibble(ENTREZID = names(gene_list), rank_stat = as.numeric(gene_list)),
      file.path(out_table_dir, paste0("ranked_gene_list_", dataset, "_", contrast, ".csv"))
    )

    run_gsea_go(gene_list, contrast, dataset, system)
  }
)

rank_summary_df <- bind_rows(rank_summary)

readr::write_csv(rank_summary_df, file.path(out_table_dir, "GSEA_ranked_gene_list_summary.csv"))
readr::write_csv(gsea_results, file.path(out_table_dir, "GO_BP_GSEA_all_contrasts.csv"))

gsea_sig <- gsea_results %>%
  filter(!is.na(p.adjust), p.adjust <= GSEA_PADJ_REPORT) %>%
  arrange(p.adjust, desc(abs(NES)))

readr::write_csv(gsea_sig, file.path(out_table_dir, "GO_BP_GSEA_significant_FDR_0.25.csv"))

# ------------------------- 5. WGCNA module ORA -------------------------------

lung_modules <- safe_read_csv(file.path(table_dir, "WGCNA_Lung_gene_modules.csv"))
tonsil_modules <- safe_read_csv(file.path(table_dir, "WGCNA_Tonsil_gene_modules.csv"))

lung_tan_genes <- lung_modules %>%
  filter(module_color == "tan" | module_numeric == 12) %>%
  pull(gene) %>%
  unique()

tonsil_blue_genes <- tonsil_modules %>%
  filter(module_color == "blue" | module_numeric == 2) %>%
  pull(gene) %>%
  unique()

module_ora_results <- bind_rows(
  run_module_ora_go(
    module_genes = lung_tan_genes,
    universe_genes = lung_modules$gene,
    module_label = "Lung tan / ME12 antiviral module",
    tissue_label = "Lung"
  ),
  run_module_ora_go(
    module_genes = tonsil_blue_genes,
    universe_genes = tonsil_modules$gene,
    module_label = "Tonsil blue / ME2 complement-coagulation module",
    tissue_label = "Tonsil"
  )
)

readr::write_csv(module_ora_results, file.path(out_table_dir, "GO_BP_WGCNA_module_ORA_all.csv"))

module_ora_sig <- module_ora_results %>%
  filter(!is.na(p.adjust), p.adjust <= ORA_PADJ_REPORT) %>%
  arrange(p.adjust)

readr::write_csv(module_ora_sig, file.path(out_table_dir, "GO_BP_WGCNA_module_ORA_significant_FDR_0.25.csv"))

module_gene_summary <- tibble(
  tissue = c("Lung", "Tonsil"),
  module = c("tan / ME12", "blue / ME2"),
  genes_in_module = c(length(lung_tan_genes), length(tonsil_blue_genes)),
  genes_mapped_to_entrez = c(
    n_distinct(map_symbols_to_entrez(lung_tan_genes)$ENTREZID),
    n_distinct(map_symbols_to_entrez(tonsil_blue_genes)$ENTREZID)
  ),
  universe_genes = c(n_distinct(lung_modules$gene), n_distinct(tonsil_modules$gene)),
  universe_mapped_to_entrez = c(
    n_distinct(map_symbols_to_entrez(lung_modules$gene)$ENTREZID),
    n_distinct(map_symbols_to_entrez(tonsil_modules$gene)$ENTREZID)
  )
)

readr::write_csv(module_gene_summary, file.path(out_table_dir, "WGCNA_module_enrichment_gene_mapping_summary.csv"))

# ------------------------------ 6. Figures -----------------------------------

if (nrow(gsea_sig) > 0) {
  gsea_plot_df <- gsea_sig %>%
    group_by(system, contrast) %>%
    slice_min(order_by = p.adjust, n = 6, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      contrast_label = str_replace_all(contrast, "_", " "),
      term_label = fct_reorder(term_label, NES),
      neg_log10_fdr = -log10(p.adjust + 1e-300)
    )

  p_gsea <- ggplot(
    gsea_plot_df,
    aes(x = NES, y = term_label, size = neg_log10_fdr, color = direction)
  ) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    facet_grid(system ~ contrast_label, scales = "free_y", space = "free_y") +
    scale_color_manual(
      values = c(
        "Enriched among upregulated genes" = "#B2182B",
        "Enriched among downregulated genes" = "#2166AC"
      )
    ) +
    scale_size_continuous(range = c(1.5, 6), name = "-log10(FDR)") +
    labs(
      title = "GO Biological Process GSEA across Nipah host-response contrasts",
      x = "Normalized enrichment score",
      y = NULL,
      color = NULL
    ) +
    theme_manuscript(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 7.5),
      strip.text.x = element_text(size = 8, angle = 0),
      strip.text.y = element_text(size = 8),
      legend.position = "bottom"
    )

  ggsave(
    file.path(out_figure_dir, "GO_BP_GSEA_summary_dotplot.png"),
    p_gsea,
    width = 18,
    height = 12,
    dpi = 400,
    bg = "white"
  )
}

if (nrow(module_ora_sig) > 0) {
  module_plot_df <- module_ora_sig %>%
    group_by(tissue, module) %>%
    slice_min(order_by = p.adjust, n = 12, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      neg_log10_fdr = -log10(p.adjust + 1e-300),
      term_label = fct_reorder(term_label, neg_log10_fdr)
    )

  p_module <- ggplot(
    module_plot_df,
    aes(x = neg_log10_fdr, y = term_label, fill = tissue)
  ) +
    geom_col(width = 0.72, color = "grey20", linewidth = 0.15) +
    facet_wrap(~ module, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Lung" = "#2C7BB6", "Tonsil" = "#D7191C")) +
    labs(
      title = "GO Biological Process enrichment of prioritized WGCNA modules",
      x = "-log10 adjusted P value",
      y = NULL,
      fill = "Tissue"
    ) +
    theme_manuscript(base_size = 10) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 10)
    )

  ggsave(
    file.path(out_figure_dir, "GO_BP_WGCNA_module_ORA_barplot.png"),
    p_module,
    width = 10,
    height = 9,
    dpi = 400,
    bg = "white"
  )
}

gsea_key_terms <- gsea_results %>%
  filter(!is.na(p.adjust)) %>%
  mutate(
    biological_theme = case_when(
      str_detect(str_to_lower(Description), "interferon|virus|viral|defense response|innate immune|cytokine") ~ "Antiviral/interferon/inflammatory",
      str_detect(str_to_lower(Description), "complement") ~ "Complement",
      str_detect(str_to_lower(Description), "coagulation|hemostasis|platelet|fibrin") ~ "Coagulation/hemostasis",
      str_detect(str_to_lower(Description), "endothelial|vascul|angiogenesis|leukocyte") ~ "Vascular/endothelial/leukocyte",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(dataset, system, contrast, biological_theme) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(dataset, contrast, biological_theme, p.adjust)

readr::write_csv(gsea_key_terms, file.path(out_table_dir, "GO_BP_GSEA_key_biological_themes.csv"))

module_key_terms <- module_ora_results %>%
  filter(!is.na(p.adjust)) %>%
  mutate(
    biological_theme = case_when(
      str_detect(str_to_lower(Description), "interferon|virus|viral|defense response|innate immune|cytokine") ~ "Antiviral/interferon/inflammatory",
      str_detect(str_to_lower(Description), "complement") ~ "Complement",
      str_detect(str_to_lower(Description), "coagulation|hemostasis|platelet|fibrin") ~ "Coagulation/hemostasis",
      str_detect(str_to_lower(Description), "endothelial|vascul|angiogenesis|leukocyte") ~ "Vascular/endothelial/leukocyte",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(tissue, module, biological_theme) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(tissue, module, biological_theme, p.adjust)

readr::write_csv(module_key_terms, file.path(out_table_dir, "GO_BP_WGCNA_module_ORA_key_biological_themes.csv"))

# ------------------------- 7. Supplemental workbook --------------------------

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb_path <- file.path(out_table_dir, "GSEA_module_enrichment_results.xlsx")
  wb <- openxlsx::createWorkbook()

  add_sheet <- function(wb, sheet, df) {
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 2)
    openxlsx::setColWidths(wb, sheet, cols = 1:min(ncol(df), 20), widths = "auto")
  }

  add_sheet(wb, "S11_GSEA_rank_summary", rank_summary_df)
  add_sheet(wb, "S12_GSEA_all", gsea_results)
  add_sheet(wb, "S13_GSEA_FDR025", gsea_sig)
  add_sheet(wb, "S14_Module_ORA_all", module_ora_results)
  add_sheet(wb, "S15_Module_ORA_FDR025", module_ora_sig)
  add_sheet(wb, "S16_Module_gene_mapping", module_gene_summary)
  add_sheet(wb, "S17_GSEA_key_themes", gsea_key_terms)
  add_sheet(wb, "S18_Module_key_themes", module_key_terms)

  openxlsx::saveWorkbook(wb, wb_path, overwrite = TRUE)
}

# ------------------------------ 8. Session -----------------------------------

version_table <- tibble(
  package = c(
    "R", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot",
    "ggplot2", "dplyr", "readr", "tidyverse", "patchwork", "openxlsx"
  ),
  version = c(
    as.character(getRversion()),
    as.character(packageVersion("clusterProfiler")),
    as.character(packageVersion("org.Hs.eg.db")),
    as.character(packageVersion("AnnotationDbi")),
    as.character(packageVersion("enrichplot")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("readr")),
    as.character(packageVersion("tidyverse")),
    as.character(packageVersion("patchwork")),
    if (requireNamespace("openxlsx", quietly = TRUE)) as.character(packageVersion("openxlsx")) else NA_character_
  )
)

readr::write_csv(version_table, file.path(out_table_dir, "GSEA_module_enrichment_package_versions.csv"))
capture.output(sessionInfo(), file = file.path(out_table_dir, "sessionInfo_GSEA_module_enrichment.txt"))

message("GSEA/module enrichment complete.")
message("Tables: ", out_table_dir)
message("Figures: ", out_figure_dir)

