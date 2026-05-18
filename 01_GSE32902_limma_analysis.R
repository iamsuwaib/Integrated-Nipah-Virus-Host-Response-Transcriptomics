# GSE32902 Nipah virus HUVEC transcriptome analysis
# Dataset: Analysis of the early transcriptome signature of Nipah virus infection
# Platform: GPL2895 CodeLink Human Whole Genome Bioarray
# Samples: MOCK_1, NiV_1, MOCK_2, NiV_2

options(stringsAsFactors = FALSE)

# ----------------------------- 1. User settings ------------------------------

BASE_DIR <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics/GSE32902"

SAMPLE_DIR <- file.path(BASE_DIR, "GSE32902_family.xml")
PLATFORM_FILE <- file.path(SAMPLE_DIR, "GPL2895-tbl-1.txt")

OUT_DIR <- file.path(BASE_DIR, "results")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ADJ_P_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1
MIN_GOOD_FLAGS <- 2

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
  "GSM813067", "NiV_2",       "NiV",      2,          "GSM813067-tbl-1.txt"
) %>%
  dplyr::mutate(
    file_path = file.path(SAMPLE_DIR, file_name),
    condition = factor(condition, levels = c("Mock", "NiV"))
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
    col.names = c(
      "ID_REF",
      paste0(sample_id, "_raw"),
      paste0(sample_id, "_value"),
      paste0(sample_id, "_flag")
    )
  ) %>%
    dplyr::mutate(ID_REF = as.character(ID_REF))
}

sample_tables <- purrr::map2(
  sample_info$file_path,
  sample_info$sample_id,
  read_sample_table
)

merged <- purrr::reduce(sample_tables, dplyr::full_join, by = "ID_REF")

expr <- merged %>%
  dplyr::select(ID_REF, dplyr::ends_with("_value"))

flag <- merged %>%
  dplyr::select(ID_REF, dplyr::ends_with("_flag"))

expr_mat <- expr %>%
  tibble::column_to_rownames("ID_REF") %>%
  as.matrix()

colnames(expr_mat) <- stringr::str_remove(colnames(expr_mat), "_value$")
expr_mat <- expr_mat[, sample_info$sample_id]
storage.mode(expr_mat) <- "numeric"

flag_mat <- flag %>%
  tibble::column_to_rownames("ID_REF") %>%
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
  dplyr::mutate(
    ID = as.character(ID),
    GB_FIRST = ifelse(
      is.na(GB_LIST) | GB_LIST == "",
      NA,
      stringr::str_split_fixed(GB_LIST, " /// |;|,", 2)[, 1]
    ),
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
  dplyr::mutate(
    SYMBOL = unname(symbol_map[GB_ACC_NOVERSION]),
    ENTREZID = unname(entrez_map[GB_ACC_NOVERSION])
  )

# ----------------------- 5. Filtering and transformation ----------------------

is_discovery <- platform$ID %in% rownames(expr_mat) &
  platform$PROBE_TYPE == "DISCOVERY"

discovery_ids <- platform$ID[is_discovery]

expr_mat <- expr_mat[rownames(expr_mat) %in% discovery_ids, ]
flag_mat <- flag_mat[rownames(flag_mat) %in% discovery_ids, ]

good_count <- rowSums(flag_mat == "G", na.rm = TRUE)
keep_quality <- good_count >= MIN_GOOD_FLAGS
keep_positive <- rowSums(expr_mat > 0, na.rm = TRUE) == ncol(expr_mat)

expr_filt <- expr_mat[keep_quality & keep_positive, ]
flag_filt <- flag_mat[keep_quality & keep_positive, ]

message("Initial probes in sample tables: ", nrow(expr))
message("Discovery probes retained: ", length(discovery_ids))
message("Probes retained after quality + positive-value filtering: ", nrow(expr_filt))

log2_expr <- log2(expr_filt + 1)

log2_expr_norm <- limma::normalizeBetweenArrays(
  log2_expr,
  method = "quantile"
)

colnames(log2_expr_norm) <- sample_info$sample_title

sample_info <- sample_info %>%
  dplyr::mutate(
    sample_title = factor(sample_title, levels = colnames(log2_expr_norm))
  )

write.csv(
  sample_info,
  file.path(TABLE_DIR, "GSE32902_sample_metadata.csv"),
  row.names = FALSE
)

# ------------------------------ 6. QC figures --------------------------------

png(
  file.path(FIG_DIR, "QC_boxplot_log2_expression.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)

boxplot(
  log2_expr_norm,
  las = 2,
  col = c("#8FB9AA", "#D95F59", "#8FB9AA", "#D95F59"),
  main = "GSE32902 normalized log2 expression",
  ylab = "log2 normalized intensity"
)

dev.off()

png(
  file.path(FIG_DIR, "QC_density_log2_expression.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)

plotDensities(
  log2_expr_norm,
  group = sample_info$condition,
  main = "Density plot: normalized log2 expression",
  legend = "topright"
)

dev.off()

pca <- prcomp(t(log2_expr_norm), scale. = TRUE)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("sample_title") %>%
  dplyr::left_join(sample_info, by = "sample_title")

percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample_title)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 4, show.legend = FALSE) +
  scale_color_manual(values = c(Mock = "#2A9D8F", NiV = "#D95F59")) +
  labs(
    title = "PCA of GSE32902 HUVEC samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(FIG_DIR, "QC_PCA.png"),
  p_pca,
  width = 6.5,
  height = 5,
  dpi = 300
)

# --------------------------- 7. Limma analysis --------------------------------

condition <- factor(sample_info$condition, levels = c("Mock", "NiV"))

design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
rownames(design) <- sample_info$sample_title

contrast_matrix <- makeContrasts(
  NiV_vs_Mock = NiV - Mock,
  levels = design
)

fit <- lmFit(log2_expr_norm, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

limma_all <- topTable(
  fit2,
  coef = "NiV_vs_Mock",
  number = Inf,
  adjust.method = "BH"
) %>%
  tibble::rownames_to_column("ID")

limma_annotated <- limma_all %>%
  dplyr::left_join(platform, by = "ID") %>%
  dplyr::relocate(
    ID,
    PROBE_NAME,
    SYMBOL,
    ENTREZID,
    GB_FIRST,
    GB_ACC_NOVERSION,
    GB_LIST,
    GI_LIST,
    PROBE_TYPE,
    .before = logFC
  ) %>%
  dplyr::arrange(adj.P.Val)

deg <- limma_annotated %>%
  dplyr::filter(adj.P.Val < ADJ_P_CUTOFF, abs(logFC) >= LOGFC_CUTOFF)

write.csv(
  limma_annotated,
  file.path(TABLE_DIR, "GSE32902_limma_all_probes_annotated.csv"),
  row.names = FALSE
)

write.csv(
  deg,
  file.path(TABLE_DIR, "GSE32902_limma_DEG_adjP_0.05_logFC_1.csv"),
  row.names = FALSE
)

write.csv(
  limma_annotated %>% dplyr::slice_head(n = 100),
  file.path(TABLE_DIR, "GSE32902_top100_ranked_genes.csv"),
  row.names = FALSE
)

message(
  "DE probes at adj.P.Val < ",
  ADJ_P_CUTOFF,
  " and |logFC| >= ",
  LOGFC_CUTOFF,
  ": ",
  nrow(deg)
)

# ----------------------------- 8. Volcano plot --------------------------------

volcano_df <- limma_annotated %>%
  dplyr::mutate(
    neg_log10_adjP = -log10(adj.P.Val),
    direction = dplyr::case_when(
      adj.P.Val < ADJ_P_CUTOFF & logFC >= LOGFC_CUTOFF ~ "Up in NiV",
      adj.P.Val < ADJ_P_CUTOFF & logFC <= -LOGFC_CUTOFF ~ "Down in NiV",
      TRUE ~ "Not significant"
    ),
    label = ifelse(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, PROBE_NAME)
  )

top_labels <- volcano_df %>%
  dplyr::filter(direction != "Not significant") %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 15)

p_volcano <- ggplot(volcano_df, aes(logFC, neg_log10_adjP, color = direction)) +
  geom_point(alpha = 0.75, size = 1.4) +
  geom_vline(
    xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_hline(
    yintercept = -log10(ADJ_P_CUTOFF),
    linetype = "dashed",
    color = "grey45"
  ) +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = label),
    size = 3.2,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Up in NiV" = "#D95F59",
      "Down in NiV" = "#2A9D8F",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    title = "GSE32902: NiV-infected HUVEC vs Mock",
    x = "log2 fold-change",
    y = "-log10 adjusted p value",
    color = ""
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(FIG_DIR, "Volcano_NiV_vs_Mock.png"),
  p_volcano,
  width = 7,
  height = 5.5,
  dpi = 300
)

# ------------------------------ 9. Heatmaps -----------------------------------

annotation_col <- sample_info %>%
  dplyr::select(sample_title, condition) %>%
  tibble::column_to_rownames("sample_title")

ann_colors <- list(
  condition = c(Mock = "#2A9D8F", NiV = "#D95F59")
)

top50_ids <- limma_annotated %>%
  dplyr::filter(!is.na(ID)) %>%
  head(50) %>%
  dplyr::pull(ID)

top50_mat <- log2_expr_norm[top50_ids, , drop = FALSE]

rownames(top50_mat) <- limma_annotated$SYMBOL[
  match(rownames(top50_mat), limma_annotated$ID)
]

rownames(top50_mat) <- ifelse(
  is.na(rownames(top50_mat)) | rownames(top50_mat) == "",
  top50_ids,
  rownames(top50_mat)
)

png(
  file.path(FIG_DIR, "Heatmap_top50_DEG.png"),
  width = 7.5,
  height = 9,
  units = "in",
  res = 300
)

pheatmap(
  top50_mat,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  fontsize_row = 7,
  main = "Top 50 ranked probes: NiV vs Mock"
)

dev.off()

antiviral_genes <- c(
  "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2",
  "OAS1", "OAS2", "OAS3", "RSAD2", "DDX58", "IFIH1",
  "IRF7", "STAT1", "STAT2", "CXCL10", "CCL5",
  "IL6", "ICAM1", "VCAM1", "SELE"
)

antiviral_ids <- limma_annotated %>%
  dplyr::filter(SYMBOL %in% antiviral_genes) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::arrange(adj.P.Val, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::pull(ID)

if (length(antiviral_ids) >= 2) {
  
  antiviral_mat <- log2_expr_norm[antiviral_ids, , drop = FALSE]
  
  rownames(antiviral_mat) <- limma_annotated$SYMBOL[
    match(rownames(antiviral_mat), limma_annotated$ID)
  ]
  
  png(
    file.path(FIG_DIR, "Heatmap_antiviral_genes.png"),
    width = 6.5,
    height = 7,
    units = "in",
    res = 300
  )
  
  pheatmap(
    antiviral_mat,
    scale = "row",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    fontsize_row = 9,
    main = "Selected antiviral and endothelial-response genes"
  )
  
  dev.off()
  
} else {
  message("Too few antiviral genes mapped for antiviral heatmap.")
}

# ---------------------- 10. Optional enrichment analysis ----------------------

if (has_clusterProfiler) {
  
  suppressPackageStartupMessages(library(clusterProfiler))
  
  entrez_deg <- deg %>%
    dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
    dplyr::pull(ENTREZID) %>%
    unique()
  
  entrez_background <- limma_annotated %>%
    dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
    dplyr::pull(ENTREZID) %>%
    unique()
  
  if (length(entrez_deg) >= 10) {
    
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
    
    write.csv(
      go_table,
      file.path(TABLE_DIR, "GSE32902_GO_BP_enrichment_DEG.csv"),
      row.names = FALSE
    )
    
    if (nrow(go_table) > 0) {
      
      png(
        file.path(FIG_DIR, "GO_BP_dotplot_DEG.png"),
        width = 8,
        height = 8,
        units = "in",
        res = 300
      )
      
      print(
        dotplot(ego_bp, showCategory = 20) +
          ggtitle("GO Biological Process enrichment: GSE32902 DEGs")
      )
      
      dev.off()
    }
    
  } else {
    message("Fewer than 10 DEG Entrez IDs mapped; skipping GO enrichment.")
  }
  
} else {
  message("clusterProfiler not available; skipping enrichment analysis.")
}

# ----------------------------- 11. Save objects -------------------------------

saveRDS(
  list(
    sample_info = sample_info,
    platform = platform,
    expression_log2_normalized = log2_expr_norm,
    flags = flag_filt,
    limma_fit = fit2,
    limma_results = limma_annotated,
    deg = deg
  ),
  file.path(OUT_DIR, "GSE32902_analysis_objects.rds")
)

sessionInfo_text <- capture.output(sessionInfo())
writeLines(sessionInfo_text, file.path(OUT_DIR, "sessionInfo.txt"))

message("Analysis complete.")
message("Results folder: ", OUT_DIR)