###############################################################################
# Integrated Nipah virus host-response signature across datasets
#
# Purpose:
#   Combine completed differential-expression outputs from:
#     - GSE32902: HUVEC, Nipah virus vs Mock
#     - GSE33133: HUVEC, NiV / NiV-dC comparisons
#     - GSE310471: African green monkey lung and tonsil, 3/4/5 DPI vs baseline
#
# Outputs:
#   tables/
#     integrated_signature_long.csv
#     integrated_signature_log2FC_matrix.csv
#     integrated_signature_padj_matrix.csv
#     integrated_signature_presence_significance.csv
#     integrated_signature_top_conserved_genes.csv
#
#   figures/
#     integrated_HUVEC_in_vivo_signature_heatmap.png
#     integrated_signature_significance_dotplot.png
#     integrated_conserved_signature_lollipop.png
#
# Notes:
#   This script assumes that individual analyses are already complete.
#   It does not run limma or DESeq2 again.
###############################################################################

options(stringsAsFactors = FALSE)

###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "tidyverse",
  "pheatmap",
  "RColorBrewer",
  "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

###############################################################################
# 2. Paths
###############################################################################

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"

gse32902_table <- file.path(
  project_dir,
  "GSE32902/results/tables/GSE32902_limma_all_probes_annotated.csv"
)

gse33133_tables <- c(
  GSE33133_HUVEC_NiV_vs_Mock =
    file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiV_vs_Mock_annotated.csv"),
  GSE33133_HUVEC_NiVdC_vs_Mock =
    file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_Mock_annotated.csv"),
  GSE33133_HUVEC_NiVdC_vs_NiV =
    file.path(project_dir, "GSE33133/results/tables/GSE33133_limma_all_NiVdC_vs_NiV_annotated.csv")
)

gse310471_tables <- c(
  GSE310471_Lung_3DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_3DPI_vs_baseline.csv"),
  GSE310471_Lung_4DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_4DPI_vs_baseline.csv"),
  GSE310471_Lung_5DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Lung_5DPI_vs_baseline.csv"),
  GSE310471_Tonsil_3DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_3DPI_vs_baseline.csv"),
  GSE310471_Tonsil_4DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_4DPI_vs_baseline.csv"),
  GSE310471_Tonsil_5DPI_vs_baseline =
    file.path(project_dir, "GSE310471/results/tables/GSE310471_DESeq2_all_Tonsil_5DPI_vs_baseline.csv")
)

out_dir <- file.path(project_dir, "integrated_results")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 3. Signature gene sets
###############################################################################

signature_categories <- tribble(
  ~gene, ~category,
  "MX1", "ISG_antiviral",
  "MX2", "ISG_antiviral",
  "OAS1", "ISG_antiviral",
  "OAS2", "ISG_antiviral",
  "OAS3", "ISG_antiviral",
  "OASL", "ISG_antiviral",
  "IFIT1", "ISG_antiviral",
  "IFIT2", "ISG_antiviral",
  "IFIT3", "ISG_antiviral",
  "IFIT5", "ISG_antiviral",
  "IFITM1", "ISG_antiviral",
  "IFIH1", "ISG_antiviral",
  "DDX58", "ISG_antiviral",
  "IRF7", "ISG_antiviral",
  "IRF9", "ISG_antiviral",
  "STAT1", "ISG_antiviral",
  "STAT2", "ISG_antiviral",
  "RSAD2", "ISG_antiviral",
  "ISG15", "ISG_antiviral",
  "HERC5", "ISG_antiviral",
  "HERC6", "ISG_antiviral",
  "SAMD9", "ISG_antiviral",
  "SAMD9L", "ISG_antiviral",
  "PARP9", "ISG_antiviral",
  "USP18", "ISG_antiviral",
  "TRIM21", "ISG_antiviral",
  "SIGLEC1", "ISG_antiviral",
  "CXCL10", "cytokine_chemokine",
  "CXCL11", "cytokine_chemokine",
  "CXCL9", "cytokine_chemokine",
  "CCL2", "cytokine_chemokine",
  "CCL5", "cytokine_chemokine",
  "IL6", "cytokine_chemokine",
  "TNF", "cytokine_chemokine",
  "ICAM1", "endothelial_activation",
  "VCAM1", "endothelial_activation",
  "SELE", "endothelial_activation",
  "ANGPT2", "endothelial_activation",
  "VWF", "endothelial_activation",
  "C3", "complement",
  "C4A", "complement",
  "C4B", "complement",
  "CFB", "complement",
  "CFD", "complement",
  "CFH", "complement",
  "CFI", "complement",
  "C1QA", "complement",
  "C1QB", "complement",
  "C1QC", "complement",
  "SERPING1", "complement",
  "F3", "coagulation",
  "FGA", "coagulation",
  "FGB", "coagulation",
  "FGG", "coagulation",
  "THBD", "coagulation",
  "PROS1", "coagulation",
  "PROC", "coagulation",
  "SERPINE1", "coagulation",
  "PLAU", "coagulation",
  "PLAUR", "coagulation"
) %>%
  distinct(gene, .keep_all = TRUE)

signature_genes <- signature_categories$gene

###############################################################################
# 4. Helper functions
###############################################################################

stop_if_missing <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      "The following required input files are missing:\n",
      paste(missing_paths, collapse = "\n")
    )
  }
}

parse_contrast_metadata <- function(contrast_name) {
  dataset <- str_extract(contrast_name, "^GSE[0-9]+")
  tissue <- case_when(
    str_detect(contrast_name, "HUVEC") ~ "HUVEC",
    str_detect(contrast_name, "Lung") ~ "Lung",
    str_detect(contrast_name, "Tonsil") ~ "Tonsil",
    TRUE ~ NA_character_
  )
  model <- case_when(
    dataset %in% c("GSE32902", "GSE33133") ~ "Human endothelial cell",
    dataset == "GSE310471" ~ "African green monkey in vivo",
    TRUE ~ NA_character_
  )
  timepoint <- case_when(
    str_detect(contrast_name, "3DPI") ~ "3 DPI",
    str_detect(contrast_name, "4DPI") ~ "4 DPI",
    str_detect(contrast_name, "5DPI") ~ "5 DPI",
    TRUE ~ "early HUVEC"
  )
  comparison <- str_remove(contrast_name, "^GSE[0-9]+_")
  tibble(
    contrast = contrast_name,
    dataset = dataset,
    model = model,
    tissue = tissue,
    timepoint = timepoint,
    comparison = comparison
  )
}

read_limma_signature <- function(path, contrast_name) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      SYMBOL = str_trim(SYMBOL),
      SYMBOL = na_if(SYMBOL, ""),
      log2FC = logFC,
      pvalue = P.Value,
      padj = adj.P.Val
    ) %>%
    filter(!is.na(SYMBOL), SYMBOL %in% signature_genes) %>%
    arrange(SYMBOL, padj, desc(abs(log2FC))) %>%
    group_by(SYMBOL) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene = SYMBOL,
      contrast = contrast_name,
      log2FC = as.numeric(log2FC),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj),
      source_id = as.character(ID)
    ) %>%
    left_join(parse_contrast_metadata(contrast_name), by = "contrast")
}

read_deseq2_signature <- function(path, contrast_name) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      Gene = str_trim(Gene),
      log2FC = log2FoldChange
    ) %>%
    filter(!is.na(Gene), Gene %in% signature_genes) %>%
    arrange(Gene, padj, desc(abs(log2FC))) %>%
    group_by(Gene) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene = Gene,
      contrast = contrast_name,
      log2FC = as.numeric(log2FC),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj),
      source_id = Gene
    ) %>%
    left_join(parse_contrast_metadata(contrast_name), by = "contrast")
}

make_named_annotation_colors <- function(annotation_df) {
  list(
    dataset = c(
      GSE32902 = "#4E79A7",
      GSE33133 = "#59A14F",
      GSE310471 = "#E15759"
    ),
    tissue = c(
      HUVEC = "#76B7B2",
      Lung = "#F28E2B",
      Tonsil = "#B07AA1"
    ),
    model = c(
      `Human endothelial cell` = "#9C755F",
      `African green monkey in vivo` = "#EDC948"
    ),
    category = c(
      ISG_antiviral = "#4E79A7",
      cytokine_chemokine = "#E15759",
      endothelial_activation = "#F28E2B",
      complement = "#59A14F",
      coagulation = "#B07AA1"
    )
  )
}

###############################################################################
# 5. Read and combine completed result tables
###############################################################################

all_input_files <- c(gse32902_table, unname(gse33133_tables), unname(gse310471_tables))
stop_if_missing(all_input_files)

gse32902_long <- read_limma_signature(
  gse32902_table,
  "GSE32902_HUVEC_NiV_vs_Mock"
)

gse33133_long <- purrr::imap_dfr(
  gse33133_tables,
  ~ read_limma_signature(path = .x, contrast_name = .y)
)

gse310471_long <- purrr::imap_dfr(
  gse310471_tables,
  ~ read_deseq2_signature(path = .x, contrast_name = .y)
)

integrated_long <- bind_rows(
  gse32902_long,
  gse33133_long,
  gse310471_long
) %>%
  left_join(signature_categories, by = "gene") %>%
  mutate(
    category = replace_na(category, "unclassified"),
    significant = !is.na(padj) & padj < 0.05 & abs(log2FC) >= 1,
    direction = case_when(
      significant & log2FC > 0 ~ "up",
      significant & log2FC < 0 ~ "down",
      TRUE ~ "not_significant"
    ),
    neg_log10_padj = case_when(
      is.na(padj) ~ NA_real_,
      padj == 0 ~ -log10(.Machine$double.xmin),
      TRUE ~ -log10(padj)
    )
  ) %>%
  arrange(category, gene, contrast)

contrast_order <- c(
  "GSE32902_HUVEC_NiV_vs_Mock",
  "GSE33133_HUVEC_NiV_vs_Mock",
  "GSE33133_HUVEC_NiVdC_vs_Mock",
  "GSE33133_HUVEC_NiVdC_vs_NiV",
  "GSE310471_Lung_3DPI_vs_baseline",
  "GSE310471_Lung_4DPI_vs_baseline",
  "GSE310471_Lung_5DPI_vs_baseline",
  "GSE310471_Tonsil_3DPI_vs_baseline",
  "GSE310471_Tonsil_4DPI_vs_baseline",
  "GSE310471_Tonsil_5DPI_vs_baseline"
)

# GSE32902 and GSE33133 are not independent HUVEC experiments: GSE33133_HUVEC_NiV_vs_Mock
# re-lists the identical four GEO sample records (GSM813064-GSM813067) already represented by
# GSE32902_HUVEC_NiV_vs_Mock (see Response to Reviewers, R2-1). GSE33133_HUVEC_NiV_vs_Mock is
# therefore retained in the raw integrated_long table, heatmap, and dotplot below for
# transparency, but is excluded from cross-dataset recurrence counting and priority-score
# calculations so that this shared evidence is not double-counted. GSE33133's NiV-dC contrasts
# are unique to GSE33133 (samples GSM814496/GSM814497) and are unaffected.
scoring_excluded_contrasts <- c("GSE33133_HUVEC_NiV_vs_Mock")

gene_order <- signature_categories %>%
  semi_join(integrated_long, by = "gene") %>%
  arrange(
    factor(
      category,
      levels = c(
        "ISG_antiviral",
        "cytokine_chemokine",
        "endothelial_activation",
        "complement",
        "coagulation"
      )
    ),
    gene
  ) %>%
  pull(gene)

integrated_long <- integrated_long %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    contrast = factor(contrast, levels = contrast_order),
    scoring_included = !as.character(contrast) %in% scoring_excluded_contrasts
  )

###############################################################################
# 6. Output integrated tables
###############################################################################

readr::write_csv(
  integrated_long,
  file.path(table_dir, "integrated_signature_long.csv")
)

log2fc_matrix <- integrated_long %>%
  select(gene, contrast, log2FC) %>%
  mutate(gene = as.character(gene), contrast = as.character(contrast)) %>%
  pivot_wider(names_from = contrast, values_from = log2FC) %>%
  arrange(match(gene, gene_order))

padj_matrix <- integrated_long %>%
  select(gene, contrast, padj) %>%
  mutate(gene = as.character(gene), contrast = as.character(contrast)) %>%
  pivot_wider(names_from = contrast, values_from = padj) %>%
  arrange(match(gene, gene_order))

presence_significance <- integrated_long %>%
  mutate(
    presence_call = case_when(
      is.na(log2FC) ~ "not_detected",
      direction == "up" ~ "significant_up",
      direction == "down" ~ "significant_down",
      TRUE ~ "detected_not_significant"
    )
  ) %>%
  select(
    gene, category, contrast, dataset, model, tissue, timepoint,
    log2FC, pvalue, padj, significant, direction, presence_call
  )

top_conserved <- integrated_long %>%
  mutate(gene = as.character(gene)) %>%
  filter(scoring_included) %>%
  group_by(gene, category) %>%
  summarise(
    n_detected = sum(!is.na(log2FC)),
    n_significant = sum(significant, na.rm = TRUE),
    n_significant_up = sum(direction == "up", na.rm = TRUE),
    n_significant_down = sum(direction == "down", na.rm = TRUE),
    n_huvec_significant_up = sum(
      tissue == "HUVEC" & direction == "up",
      na.rm = TRUE
    ),
    n_in_vivo_significant_up = sum(
      tissue %in% c("Lung", "Tonsil") & direction == "up",
      na.rm = TRUE
    ),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    min_padj = ifelse(is.infinite(min_padj), NA_real_, min_padj),
    conserved_call = case_when(
      n_huvec_significant_up >= 1 & n_in_vivo_significant_up >= 2 ~
        "conserved_HUVEC_and_in_vivo",
      n_huvec_significant_up >= 1 & n_in_vivo_significant_up < 2 ~
        "mainly_HUVEC",
      n_huvec_significant_up == 0 & n_in_vivo_significant_up >= 2 ~
        "mainly_in_vivo",
      TRUE ~ "not_strongly_conserved"
    )
  ) %>%
  arrange(
    desc(n_significant_up),
    desc(n_huvec_significant_up),
    desc(n_in_vivo_significant_up),
    min_padj
  )

readr::write_csv(
  log2fc_matrix,
  file.path(table_dir, "integrated_signature_log2FC_matrix.csv")
)

readr::write_csv(
  padj_matrix,
  file.path(table_dir, "integrated_signature_padj_matrix.csv")
)

readr::write_csv(
  presence_significance,
  file.path(table_dir, "integrated_signature_presence_significance.csv")
)

readr::write_csv(
  top_conserved,
  file.path(table_dir, "integrated_signature_top_conserved_genes.csv")
)

###############################################################################
# 7. Integrated heatmap
###############################################################################

heatmap_mat <- log2fc_matrix %>%
  column_to_rownames("gene") %>%
  as.matrix()

heatmap_mat <- heatmap_mat[, contrast_order[contrast_order %in% colnames(heatmap_mat)], drop = FALSE]
heatmap_mat <- heatmap_mat[gene_order[gene_order %in% rownames(heatmap_mat)], , drop = FALSE]

heatmap_mat_capped <- pmax(pmin(heatmap_mat, 5), -5)

row_annotation <- signature_categories %>%
  filter(gene %in% rownames(heatmap_mat_capped)) %>%
  arrange(match(gene, rownames(heatmap_mat_capped))) %>%
  column_to_rownames("gene") %>%
  select(category)

column_annotation <- map_dfr(
  colnames(heatmap_mat_capped),
  parse_contrast_metadata
) %>%
  column_to_rownames("contrast") %>%
  select(dataset, tissue, model)

annotation_colors <- make_named_annotation_colors(column_annotation)

heatmap_palette <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(11, "RdBu"))
)(101)

# Mark the GSE33133 NiV-vs-Mock column as a shared/duplicate sample set (see R2-1): it is
# displayed here for transparency but excluded from recurrence/priority-score counting.
display_labels_col <- ifelse(
  colnames(heatmap_mat_capped) %in% scoring_excluded_contrasts,
  paste0(colnames(heatmap_mat_capped), "*"),
  colnames(heatmap_mat_capped)
)

png(
  filename = file.path(figure_dir, "integrated_HUVEC_in_vivo_signature_heatmap.png"),
  width = 6800,
  height = 6200,
  res = 300
)

pheatmap::pheatmap(
  heatmap_mat_capped,
  color = heatmap_palette,
  breaks = seq(-5, 5, length.out = 102),
  border_color = NA,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annotation,
  annotation_col = column_annotation,
  annotation_colors = annotation_colors,
  labels_col = display_labels_col,
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 45,
  cellwidth = 125,
  cellheight = 22,
  legend = TRUE
)

dev.off()

###############################################################################
# 8. Significance dotplot
###############################################################################

dotplot_df <- integrated_long %>%
  mutate(
    gene = factor(as.character(gene), levels = rev(gene_order)),
    contrast = factor(as.character(contrast), levels = contrast_order),
    neg_log10_padj_plot = pmin(replace_na(neg_log10_padj, 0), 20),
    significant_label = ifelse(significant, "FDR < 0.05 and |log2FC| >= 1", "not significant")
  )

# Same duplicate-column labeling as the heatmap above (see R2-1).
dotplot_contrast_labels <- setNames(
  ifelse(
    contrast_order %in% scoring_excluded_contrasts,
    paste0(contrast_order, "*"),
    contrast_order
  ),
  contrast_order
)

png(
  filename = file.path(figure_dir, "integrated_signature_significance_dotplot.png"),
  width = 6200,
  height = 5600,
  res = 300
)

print(
  ggplot(dotplot_df, aes(x = contrast, y = gene)) +
    geom_point(
      aes(size = neg_log10_padj_plot, fill = log2FC),
      shape = 21,
      color = "grey20",
      stroke = 0.2,
      alpha = 0.9
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-5, 5),
      oob = scales::squish,
      name = "log2FC"
    ) +
    scale_size_continuous(
      range = c(1, 6),
      name = "-log10(FDR)",
      breaks = c(1.3, 2, 5, 10)
    ) +
    scale_x_discrete(labels = dotplot_contrast_labels) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 11, face = "plain"),
      strip.text.y = element_text(angle = 0, face = "bold", size = 12),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(20, 80, 80, 40)
    )
)

dev.off()

###############################################################################
# 9. Conserved signature lollipop/barplot - publication style
###############################################################################

barplot_df <- top_conserved %>%
  dplyr::filter(n_significant_up > 0) %>%
  dplyr::slice_max(order_by = n_significant_up, n = 30, with_ties = FALSE) %>%
  dplyr::arrange(n_significant_up, mean_log2FC) %>%
  dplyr::mutate(
    gene = factor(gene, levels = gene),
    conserved_call = factor(
      conserved_call,
      levels = c(
        "conserved_HUVEC_and_in_vivo",
        "mainly_in_vivo",
        "mainly_HUVEC",
        "not_strongly_conserved"
      )
    )
  )

png(
  filename = file.path(figure_dir, "integrated_conserved_signature_lollipop.png"),
  width = 4200,
  height = 3200,
  res = 300
)

print(
  ggplot(barplot_df, aes(y = gene, x = n_significant_up)) +
    geom_segment(
      aes(x = 0, xend = n_significant_up, y = gene, yend = gene),
      linewidth = 0.9,
      color = "grey70"
    ) +
    geom_point(
      aes(
        fill = conserved_call,
        size = n_in_vivo_significant_up
      ),
      shape = 21,
      color = "grey20",
      stroke = 0.4,
      alpha = 0.95
    ) +
    geom_text(
      aes(
        label = paste0(
          "HUVEC: ", n_huvec_significant_up,
          " | in vivo: ", n_in_vivo_significant_up
        )
      ),
      hjust = -0.15,
      size = 3.2,
      color = "grey20"
    ) +
    scale_fill_manual(
      values = c(
        conserved_HUVEC_and_in_vivo = "#762A83",
        mainly_in_vivo = "#E66101",
        mainly_HUVEC = "#1F78B4",
        not_strongly_conserved = "#999999"
      ),
      labels = c(
        conserved_HUVEC_and_in_vivo = "Conserved in HUVEC and in vivo",
        mainly_in_vivo = "Mainly in vivo",
        mainly_HUVEC = "Mainly HUVEC",
        not_strongly_conserved = "Not strongly conserved"
      ),
      name = NULL
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name = "In vivo significant\nupregulated contrasts"
    ) +
    scale_x_continuous(
      limits = c(0, max(barplot_df$n_significant_up) + 3),
      breaks = seq(0, max(barplot_df$n_significant_up) + 1, by = 1),
      expand = expansion(mult = c(0.01, 0.08))
    ) +
    labs(
      x = "Number of significant upregulated contrasts",
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_text(face = "plain", size = 11),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 11),
      legend.title = element_text(face = "bold", size = 11),
      plot.margin = margin(25, 90, 35, 35)
    )
)

dev.off()
###############################################################################
# 10. Console summary
###############################################################################

message("Integrated signature analysis complete.")
message("Tables written to: ", table_dir)
message("Figures written to: ", figure_dir)
message("")
message("Top conserved genes:")
print(
  top_conserved %>%
    select(
      gene,
      category,
      n_significant_up,
      n_huvec_significant_up,
      n_in_vivo_significant_up,
      conserved_call,
      min_padj
    ) %>%
    head(20)
)

