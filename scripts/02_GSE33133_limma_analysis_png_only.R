# GSE33133 Nipah virus C-protein deletion HUVEC transcriptome analysis
# Dataset: Early transcriptome signature of primary endothelial cells infected
#          by wild-type Nipah virus versus Nipah virus deleted for C protein
# Platform: GPL2895 CodeLink Human Whole Genome Bioarray
#
# Samples:
#   MOCK_1     GSM813064
#   NiV_1      GSM813065
#   MOCK_2     GSM813066
#   NiV_2      GSM813067
#   NiV_dC_1   GSM814496
#   NiV_dC_2   GSM814497
#
# Main contrasts:
#   1. NiV_vs_Mock       = wild-type Nipah virus infection vs mock
#   2. NiVdC_vs_Mock     = C-deleted Nipah virus infection vs mock
#   3. NiVdC_vs_NiV      = C-deleted Nipah virus vs wild-type Nipah virus
#
# Output style:
#   - PNG figures only. No PDFs are produced.
#
# How to use:
#   1. Open this script in RStudio.
#   2. Change BASE_DIR if needed.
#   3. Run the whole script.

options(stringsAsFactors = FALSE)

# ----------------------------- 1. User settings ------------------------------

BASE_DIR <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics/GSE33133"

SAMPLE_DIR <- file.path(BASE_DIR, "GSE33133_family.xml")
PLATFORM_FILE <- file.path(SAMPLE_DIR, "GPL2895-tbl-1.txt")

OUT_DIR <- file.path(BASE_DIR, "results")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ADJ_P_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1
MIN_GOOD_FLAGS <- 3

# ------------------------- 2. Package management -----------------------------

cran_pkgs <- c("tidyverse", "pheatmap", "ggrepel")
bioc_pkgs <- c("limma", "Biobase", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler")

install_if_missing <- function(pkgs, source = c("cran", "bioc")) {
  source <- match.arg(source)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))

  if (source == "cran") {
    install.packages(missing, dependencies = TRUE)
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
}

install_if_missing(cran_pkgs, "cran")
install_if_missing(bioc_pkgs, "bioc")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)

# --------------------------- 3. Read local files ------------------------------

stopifnot(dir.exists(BASE_DIR))
stopifnot(dir.exists(SAMPLE_DIR))
stopifnot(file.exists(PLATFORM_FILE))

sample_info <- tibble::tribble(
  ~sample_id,   ~sample_title, ~condition, ~replicate, ~file_name,
  "GSM813064", "MOCK_1",      "Mock",     1,          "GSM813064-tbl-1.txt",
  "GSM813065", "NiV_1",       "NiV",      1,          "GSM813065-tbl-1.txt",
  "GSM813066", "MOCK_2",      "Mock",     2,          "GSM813066-tbl-1.txt",
  "GSM813067", "NiV_2",       "NiV",      2,          "GSM813067-tbl-1.txt",
  "GSM814496", "NiV_dC_1",    "NiV_dC",   1,          "GSM814496-tbl-1.txt",
  "GSM814497", "NiV_dC_2",    "NiV_dC",   2,          "GSM814497-tbl-1.txt"
) %>%
  mutate(
    file_path = file.path(SAMPLE_DIR, file_name),
    condition = factor(condition, levels = c("Mock", "NiV", "NiV_dC"))
  )

if (!all(file.exists(sample_info$file_path))) {
  stop("One or more sample table files are missing. Check SAMPLE_DIR and file names.")
}

read_sample_table <- function(path, sample_id) {
  read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    col.names = c("ID_REF", paste0(sample_id, "_raw"),
                  paste0(sample_id, "_value"),
                  paste0(sample_id, "_flag"))
  ) %>%
    mutate(ID_REF = as.character(ID_REF))
}

sample_tables <- purrr::map2(sample_info$file_path, sample_info$sample_id, read_sample_table)
merged <- purrr::reduce(sample_tables, full_join, by = "ID_REF")

expr <- merged %>% dplyr::select(ID_REF, dplyr::ends_with("_value"))
flag <- merged %>% dplyr::select(ID_REF, dplyr::ends_with("_flag"))

expr_mat <- expr %>%
  column_to_rownames("ID_REF") %>%
  as.matrix()
colnames(expr_mat) <- stringr::str_remove(colnames(expr_mat), "_value$")
expr_mat <- expr_mat[, sample_info$sample_id]
storage.mode(expr_mat) <- "numeric"

flag_mat <- flag %>%
  column_to_rownames("ID_REF") %>%
  as.matrix()
colnames(flag_mat) <- stringr::str_remove(colnames(flag_mat), "_flag$")
flag_mat <- flag_mat[, sample_info$sample_id]

# ------------------------- 4. Read platform annotation ------------------------

platform <- read.delim(
  PLATFORM_FILE,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c(
    "ID", "LOGICAL_ROW", "LOGICAL_COL", "PROBE_NAME", "PROBE_TYPE",
    "PUB_PROBE_TARGETS", "SPOT_ID", "GB_LIST", "GI_LIST"
  ),
  fill = TRUE
) %>%
  mutate(
    ID = as.character(ID),
    GB_FIRST = ifelse(is.na(GB_LIST) | GB_LIST == "", NA,
                      str_split_fixed(GB_LIST, " /// |;|,", 2)[, 1]),
    GB_FIRST = stringr::str_trim(GB_FIRST),
    GB_ACC_NOVERSION = stringr::str_replace(GB_FIRST, "\\.[0-9]+$", "")
  )

acc_keys <- platform$GB_ACC_NOVERSION[
  !is.na(platform$GB_ACC_NOVERSION) & platform$GB_ACC_NOVERSION != ""
]

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(acc_keys),
  column = "SYMBOL",
  keytype = "ACCNUM",
  multiVals = "first"
)

entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(acc_keys),
  column = "ENTREZID",
  keytype = "ACCNUM",
  multiVals = "first"
)

platform <- platform %>%
  mutate(
    SYMBOL = unname(symbol_map[GB_ACC_NOVERSION]),
    ENTREZID = unname(entrez_map[GB_ACC_NOVERSION])
  )

# ----------------------- 5. Filtering and transformation ----------------------

discovery_ids <- platform %>%
  filter(PROBE_TYPE == "DISCOVERY", ID %in% rownames(expr_mat)) %>%
  pull(ID)

expr_mat <- expr_mat[rownames(expr_mat) %in% discovery_ids, ]
flag_mat <- flag_mat[rownames(flag_mat) %in% discovery_ids, ]

good_count <- rowSums(flag_mat == "G", na.rm = TRUE)
keep_quality <- good_count >= MIN_GOOD_FLAGS
keep_positive <- rowSums(expr_mat > 0, na.rm = TRUE) == ncol(expr_mat)

expr_filt <- expr_mat[keep_quality & keep_positive, ]
flag_filt <- flag_mat[keep_quality & keep_positive, ]

message("Initial probes in sample tables: ", nrow(expr))
message("Discovery probes retained before quality filtering: ", length(discovery_ids))
message("Probes retained after quality + positive-value filtering: ", nrow(expr_filt))

log2_expr <- log2(expr_filt + 1)
log2_expr_norm <- limma::normalizeBetweenArrays(log2_expr, method = "quantile")
colnames(log2_expr_norm) <- sample_info$sample_title

write.csv(sample_info, file.path(TABLE_DIR, "GSE33133_sample_metadata.csv"), row.names = FALSE)

filter_summary <- tibble(
  total_rows_in_tables = nrow(expr),
  discovery_probes = length(discovery_ids),
  min_good_flags = MIN_GOOD_FLAGS,
  probes_passing_quality = sum(keep_quality),
  probes_positive_all_samples = sum(keep_positive),
  probes_used_for_limma = nrow(expr_filt)
)
write.csv(filter_summary, file.path(TABLE_DIR, "GSE33133_filter_summary.csv"), row.names = FALSE)

# ------------------------------ 6. QC figures --------------------------------

condition_colors <- c(Mock = "#2A9D8F", NiV = "#D95F59", NiV_dC = "#6A7FDB")
sample_colors <- condition_colors[as.character(sample_info$condition)]

png(file.path(FIG_DIR, "QC_boxplot_log2_expression.png"), width = 1800, height = 1200, res = 180)
boxplot(
  log2_expr_norm,
  las = 2,
  col = sample_colors,
  main = "GSE33133 normalized log2 expression",
  ylab = "log2 normalized intensity"
)
legend("topright", legend = names(condition_colors), fill = condition_colors, bty = "n")
dev.off()

png(file.path(FIG_DIR, "QC_density_log2_expression.png"), width = 1800, height = 1200, res = 180)
plotDensities(
  log2_expr_norm,
  group = sample_info$condition,
  main = "Density plot: normalized log2 expression",
  legend = "topright"
)
dev.off()

pca <- prcomp(t(log2_expr_norm), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample_title") %>%
  left_join(sample_info, by = "sample_title")
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample_title)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 4, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA of GSE33133 HUVEC samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(FIG_DIR, "QC_PCA.png"), p_pca, width = 8, height = 6, dpi = 300)

# --------------------------- 7. Limma analysis --------------------------------

condition <- factor(sample_info$condition, levels = c("Mock", "NiV", "NiV_dC"))
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
rownames(design) <- sample_info$sample_title

contrast_matrix <- makeContrasts(
  NiV_vs_Mock = NiV - Mock,
  NiVdC_vs_Mock = NiV_dC - Mock,
  NiVdC_vs_NiV = NiV_dC - NiV,
  levels = design
)

fit <- lmFit(log2_expr_norm, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

annotate_results <- function(tt) {
  tt %>%
    rownames_to_column("ID") %>%
    left_join(platform, by = "ID") %>%
    relocate(ID, PROBE_NAME, SYMBOL, ENTREZID, GB_FIRST, GB_ACC_NOVERSION,
             GB_LIST, GI_LIST, PROBE_TYPE, .before = logFC) %>%
    arrange(adj.P.Val)
}

contrast_names <- colnames(contrast_matrix)

all_results <- list()
deg_results <- list()

for (contrast in contrast_names) {
  res <- topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH") %>%
    annotate_results()

  deg <- res %>%
    filter(adj.P.Val < ADJ_P_CUTOFF, abs(logFC) >= LOGFC_CUTOFF)

  all_results[[contrast]] <- res
  deg_results[[contrast]] <- deg

  write.csv(res,
            file.path(TABLE_DIR, paste0("GSE33133_limma_all_", contrast, "_annotated.csv")),
            row.names = FALSE)
  write.csv(deg,
            file.path(TABLE_DIR, paste0("GSE33133_DEG_", contrast,
                                        "_adjP_0.05_logFC_1.csv")),
            row.names = FALSE)
  write.csv(res %>% slice_head(n = 100),
            file.path(TABLE_DIR, paste0("GSE33133_top100_", contrast, ".csv")),
            row.names = FALSE)

  message(contrast, ": ", nrow(deg), " DE probes at adj.P.Val < ",
          ADJ_P_CUTOFF, " and |logFC| >= ", LOGFC_CUTOFF)
}

deg_summary <- tibble(
  contrast = contrast_names,
  n_deg = vapply(deg_results, nrow, integer(1)),
  n_up = vapply(deg_results, function(x) sum(x$logFC >= LOGFC_CUTOFF), integer(1)),
  n_down = vapply(deg_results, function(x) sum(x$logFC <= -LOGFC_CUTOFF), integer(1))
)
write.csv(deg_summary, file.path(TABLE_DIR, "GSE33133_DEG_summary.csv"), row.names = FALSE)

# ----------------------------- 8. Volcano plots -------------------------------

plot_volcano <- function(res, contrast) {
  volcano_df <- res %>%
    mutate(
      neg_log10_adjP = -log10(adj.P.Val),
      direction = case_when(
        adj.P.Val < ADJ_P_CUTOFF & logFC >= LOGFC_CUTOFF ~ "Up",
        adj.P.Val < ADJ_P_CUTOFF & logFC <= -LOGFC_CUTOFF ~ "Down",
        TRUE ~ "Not significant"
      ),
      label = ifelse(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, PROBE_NAME)
    )

  top_labels <- volcano_df %>%
    filter(direction != "Not significant") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 15)

  p <- ggplot(volcano_df, aes(logFC, neg_log10_adjP, color = direction)) +
    geom_point(alpha = 0.75, size = 1.4) +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
               linetype = "dashed", color = "grey45") +
    geom_hline(yintercept = -log10(ADJ_P_CUTOFF),
               linetype = "dashed", color = "grey45") +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = label),
      size = 3.2,
      max.overlaps = 30,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "Up" = "#D95F59",
      "Down" = "#2A9D8F",
      "Not significant" = "grey70"
    )) +
    labs(
      title = paste0("GSE33133: ", contrast),
      x = "log2 fold-change",
      y = "-log10 adjusted p value",
      color = ""
    ) +
    theme_bw(base_size = 13)

  ggsave(file.path(FIG_DIR, paste0("Volcano_", contrast, ".png")),
         p, width = 8, height = 6, dpi = 300)
}

purrr::walk(contrast_names, ~ plot_top_heatmap(all_results[[.x]], .x, n = 50))

# ------------------------------ 9. Heatmaps -----------------------------------

annotation_col <- sample_info %>%
  dplyr::select(sample_title, condition) %>%
  column_to_rownames("sample_title")
ann_colors <- list(condition = condition_colors)

plot_top_heatmap <- function(res, contrast, top_n = 50) {
  top_ids <- res %>%
    dplyr::filter(!is.na(ID)) %>%
    dplyr::slice_head(n = min(top_n, nrow(res))) %>%
    dplyr::pull(ID)
  
  mat <- log2_expr_norm[top_ids, , drop = FALSE]
  
  labels <- res$SYMBOL[match(rownames(mat), res$ID)]
  rownames(mat) <- ifelse(is.na(labels) | labels == "", rownames(mat), labels)
  
  png(
    file.path(FIG_DIR, paste0("Heatmap_top", top_n, "_", contrast, ".png")),
    width = 2200,
    height = 2800,
    res = 220
  )
  
  pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    fontsize_row = 7,
    main = paste0("Top ", top_n, " ranked probes: ", contrast)
  )
  
  dev.off()
}

purrr::walk(contrast_names, ~ plot_top_heatmap(all_results[[.x]], .x, top_n = 50))

antiviral_genes <- c(
  "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "ISG15", "MX1", "MX2",
  "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "DDX58", "IFIH1", "IRF7",
  "STAT1", "STAT2", "CXCL10", "CXCL11", "CCL5", "IL6", "ICAM1", "VCAM1",
  "SELE", "HERC5", "HERC6", "SAMD9", "SAMD9L", "PARP9", "USP18", "TRIM21"
)

best_by_gene <- all_results[["NiVdC_vs_NiV"]] %>%
  filter(SYMBOL %in% antiviral_genes) %>%
  group_by(SYMBOL) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

if (nrow(best_by_gene) >= 2) {
  antiviral_mat <- log2_expr_norm[best_by_gene$ID, , drop = FALSE]
  rownames(antiviral_mat) <- best_by_gene$SYMBOL

  png(file.path(FIG_DIR, "Heatmap_antiviral_endothelial_genes.png"),
      width = 1800, height = 2400, res = 220)
  pheatmap(
    antiviral_mat,
    scale = "row",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = "Selected antiviral and endothelial-response genes"
  )
  dev.off()
}

# ---------------------- 10. Contrast overlap summaries ------------------------

deg_gene_sets <- purrr::map(deg_results, function(x) {
  x %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    pull(SYMBOL) %>%
    unique()
})

overlap_summary <- expand.grid(
  contrast_1 = names(deg_gene_sets),
  contrast_2 = names(deg_gene_sets),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    n_overlap = length(intersect(deg_gene_sets[[contrast_1]],
                                 deg_gene_sets[[contrast_2]]))
  ) %>%
  ungroup()

write.csv(overlap_summary,
          file.path(TABLE_DIR, "GSE33133_DEG_symbol_overlap_summary.csv"),
          row.names = FALSE)

overlap_mat <- overlap_summary %>%
  tidyr::pivot_wider(names_from = contrast_2, values_from = n_overlap) %>%
  column_to_rownames("contrast_1") %>%
  as.matrix()

png(file.path(FIG_DIR, "Heatmap_DEG_symbol_overlap.png"), width = 1400, height = 1200, res = 180)
pheatmap(
  overlap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_color = "black",
  main = "Overlap of DEG gene symbols across contrasts"
)
dev.off()

# ---------------------- 11. Optional GO enrichment ----------------------------

run_go <- function(deg, res, contrast) {
  if (!has_clusterProfiler) return(NULL)
  suppressPackageStartupMessages(library(clusterProfiler))

  entrez_deg <- deg %>%
    filter(!is.na(ENTREZID), ENTREZID != "") %>%
    pull(ENTREZID) %>%
    unique()

  entrez_background <- res %>%
    filter(!is.na(ENTREZID), ENTREZID != "") %>%
    pull(ENTREZID) %>%
    unique()

  if (length(entrez_deg) < 10) {
    message(contrast, ": fewer than 10 mapped DEG Entrez IDs; skipping GO.")
    return(NULL)
  }

  ego_bp <- enrichGO(
    gene = entrez_deg,
    universe = entrez_background,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  go_table <- as.data.frame(ego_bp)
  write.csv(go_table,
            file.path(TABLE_DIR, paste0("GSE33133_GO_BP_", contrast, ".csv")),
            row.names = FALSE)

  if (nrow(go_table) > 0) {
    p_go <- dotplot(ego_bp, showCategory = min(20, nrow(go_table))) +
      ggtitle(paste0("GO Biological Process: ", contrast)) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold")
      )

    ggsave(file.path(FIG_DIR, paste0("GO_BP_dotplot_", contrast, ".png")),
           p_go, width = 9, height = 7, dpi = 300)
  }

  return(ego_bp)
}

go_results <- purrr::map(contrast_names, ~ run_go(deg_results[[.x]], all_results[[.x]], .x))
names(go_results) <- contrast_names

# ----------------------------- 12. Save objects -------------------------------

saveRDS(
  list(
    sample_info = sample_info,
    filter_summary = filter_summary,
    platform = platform,
    expression_log2_normalized = log2_expr_norm,
    flags = flag_filt,
    limma_fit = fit2,
    all_results = all_results,
    deg_results = deg_results,
    deg_summary = deg_summary,
    go_results = go_results
  ),
  file.path(OUT_DIR, "GSE33133_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

message("Analysis complete.")
message("Results folder: ", OUT_DIR)
