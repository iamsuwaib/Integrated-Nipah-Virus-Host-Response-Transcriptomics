# GSE310471 Nipah virus infection progression in African green monkeys
# PNG-only output

options(stringsAsFactors = FALSE)

# ----------------------------- 1. User settings ------------------------------

BASE_DIR <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics/GSE310471"

COUNTS_FILE_CANDIDATES <- c(
  file.path(BASE_DIR, "GSE310471_counts.tsv", "GSE310471_counts.tsv"),
  file.path(BASE_DIR, "GSE310471_counts.tsv.gz")
)

COUNTS_FILE <- COUNTS_FILE_CANDIDATES[file.exists(COUNTS_FILE_CANDIDATES)][1]
if (is.na(COUNTS_FILE)) {
  stop("Could not find GSE310471_counts.tsv or GSE310471_counts.tsv.gz under BASE_DIR.")
}

OUT_DIR <- file.path(BASE_DIR, "results")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
LOG2FC_CUTOFF <- 1

MIN_COUNT <- 10
MIN_SAMPLES <- 4

# ------------------------- 2. Package management -----------------------------

cran_pkgs <- c("tidyverse", "pheatmap", "ggrepel", "RColorBrewer")
bioc_pkgs <- c("DESeq2", "Biobase", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler")

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
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)

# ------------------------- 3. Build sample metadata ---------------------------

sample_ids <- sprintf("sample%03d", 59:93)

sample_info <- tibble::tibble(
  sample_id = sample_ids,
  geo_accession = paste0("GSM", 9301871:9301905),
  tissue = c(rep("Lung", 17), rep("Tonsil", 18)),
  dpi = c(
    rep("baseline", 5), rep("3DPI", 4), rep("4DPI", 4), rep("5DPI", 4),
    rep("baseline", 6), rep("3DPI", 4), rep("4DPI", 4), rep("5DPI", 4)
  ),
  replicate = c(
    LETTERS[1:5], LETTERS[1:4], LETTERS[1:4], LETTERS[1:4],
    LETTERS[1:6], LETTERS[1:4], LETTERS[1:4], LETTERS[1:4]
  )
) %>%
  dplyr::mutate(
    tissue = factor(tissue, levels = c("Lung", "Tonsil")),
    dpi = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI")),
    group = factor(paste(tissue, dpi, sep = "_")),
    sample_label = paste(tissue, dpi, replicate, sep = "_")
  )

write.csv(
  sample_info,
  file.path(TABLE_DIR, "GSE310471_sample_metadata.csv"),
  row.names = FALSE
)

# ----------------------------- 4. Read counts ---------------------------------

counts_df <- read.delim(
  COUNTS_FILE,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE
)

if (!"Gene" %in% colnames(counts_df)) {
  stop("Expected first column named 'Gene'. Please check the count file format.")
}

missing_samples <- setdiff(sample_info$sample_id, colnames(counts_df))

if (length(missing_samples) > 0) {
  stop("Count file is missing expected sample columns: ",
       paste(missing_samples, collapse = ", "))
}

counts_df <- counts_df %>%
  dplyr::select(Gene, dplyr::all_of(sample_info$sample_id)) %>%
  dplyr::mutate(Gene = make.unique(as.character(Gene)))

count_mat <- counts_df %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

storage.mode(count_mat) <- "integer"

stopifnot(all(colnames(count_mat) == sample_info$sample_id))

write.csv(
  tibble::tibble(
    n_genes_total = nrow(count_mat),
    n_samples_total = ncol(count_mat),
    count_file = COUNTS_FILE
  ),
  file.path(TABLE_DIR, "GSE310471_input_summary.csv"),
  row.names = FALSE
)

# ------------------------------- 5. QC plots ----------------------------------

library_sizes <- colSums(count_mat)

lib_df <- sample_info %>%
  dplyr::mutate(library_size = library_sizes[sample_id])

p_lib <- ggplot(lib_df, aes(sample_label, library_size / 1e6, fill = tissue)) +
  geom_col(width = 0.75) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(Lung = "#2A9D8F", Tonsil = "#6A7FDB")) +
  labs(
    title = "GSE310471 library sizes",
    x = NULL,
    y = "Library size, millions of counts"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  )

ggsave(
  file.path(FIG_DIR, "QC_library_sizes.png"),
  p_lib,
  width = 12,
  height = 5.5,
  dpi = 300
)

dds_all <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = as.data.frame(sample_info %>% tibble::column_to_rownames("sample_id")),
  design = ~ tissue + dpi
)

keep_all <- rowSums(counts(dds_all) >= MIN_COUNT) >= MIN_SAMPLES
dds_all <- dds_all[keep_all, ]

dds_all <- estimateSizeFactors(dds_all)

vst_all <- vst(dds_all, blind = TRUE)
vst_mat_all <- assay(vst_all)

write.csv(
  tibble::tibble(
    filtering_scope = "all_samples_QC",
    min_count = MIN_COUNT,
    min_samples = MIN_SAMPLES,
    genes_retained = nrow(vst_mat_all)
  ),
  file.path(TABLE_DIR, "GSE310471_QC_filter_summary.csv"),
  row.names = FALSE
)

pca_data <- plotPCA(vst_all, intgroup = c("tissue", "dpi"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_data <- pca_data %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, PC1, PC2) %>%
  dplyr::left_join(sample_info, by = "sample_id")

p_pca_all <- ggplot(
  pca_data,
  aes(PC1, PC2, color = dpi, shape = tissue, label = sample_label)
) +
  geom_point(size = 3.6) +
  ggrepel::geom_text_repel(size = 2.8, max.overlaps = 60, show.legend = FALSE) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "GSE310471 PCA: all lung and tonsil samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(FIG_DIR, "QC_PCA_all_samples.png"),
  p_pca_all,
  width = 9,
  height = 7,
  dpi = 300
)

sample_dists <- dist(t(vst_mat_all))
sample_dist_mat <- as.matrix(sample_dists)

rownames(sample_dist_mat) <- sample_info$sample_label[
  match(rownames(sample_dist_mat), sample_info$sample_id)
]

colnames(sample_dist_mat) <- sample_info$sample_label[
  match(colnames(sample_dist_mat), sample_info$sample_id)
]

annotation_col <- sample_info %>%
  dplyr::select(sample_id, tissue, dpi) %>%
  tibble::column_to_rownames("sample_id")

rownames(annotation_col) <- sample_info$sample_label

ann_colors <- list(
  tissue = c(Lung = "#2A9D8F", Tonsil = "#6A7FDB"),
  dpi = c(
    baseline = "#666666",
    `3DPI` = "#D95F59",
    `4DPI` = "#E9C46A",
    `5DPI` = "#F4A261"
  )
)

png(
  file.path(FIG_DIR, "QC_sample_distance_heatmap_all.png"),
  width = 2600,
  height = 2400,
  res = 220
)

pheatmap(
  sample_dist_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  fontsize = 7,
  main = "Sample distance heatmap: all samples"
)

dev.off()

# ------------------------- 6. DESeq2 by tissue --------------------------------

run_tissue_deseq <- function(tissue_name) {

  message("Running DESeq2 for tissue: ", tissue_name)

  tissue_samples <- sample_info %>%
    dplyr::filter(tissue == tissue_name) %>%
    droplevels()

  tissue_counts <- count_mat[, tissue_samples$sample_id, drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = tissue_counts,
    colData = as.data.frame(tissue_samples %>% tibble::column_to_rownames("sample_id")),
    design = ~ dpi
  )

  keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  vst_obj <- vst(dds, blind = FALSE)
  vst_mat <- assay(vst_obj)

  pca_tissue <- plotPCA(vst_obj, intgroup = "dpi", returnData = TRUE)
  percent_var_tissue <- round(100 * attr(pca_tissue, "percentVar"))

  pca_tissue <- pca_tissue %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::select(sample_id, PC1, PC2) %>%
    dplyr::left_join(tissue_samples, by = "sample_id")

  p <- ggplot(
    pca_tissue,
    aes(PC1, PC2, color = dpi, label = sample_label)
  ) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3.2, max.overlaps = 40, show.legend = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = paste0("GSE310471 PCA: ", tissue_name),
      x = paste0("PC1 (", percent_var_tissue[1], "%)"),
      y = paste0("PC2 (", percent_var_tissue[2], "%)")
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(FIG_DIR, paste0("QC_PCA_", tissue_name, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )

  tissue_dists <- dist(t(vst_mat))
  tissue_dist_mat <- as.matrix(tissue_dists)

  rownames(tissue_dist_mat) <- tissue_samples$sample_label[
    match(rownames(tissue_dist_mat), tissue_samples$sample_id)
  ]

  colnames(tissue_dist_mat) <- tissue_samples$sample_label[
    match(colnames(tissue_dist_mat), tissue_samples$sample_id)
  ]

  tissue_ann <- tissue_samples %>%
    dplyr::select(sample_id, dpi) %>%
    tibble::column_to_rownames("sample_id")

  rownames(tissue_ann) <- tissue_samples$sample_label

  png(
    file.path(FIG_DIR, paste0("QC_sample_distance_heatmap_", tissue_name, ".png")),
    width = 2000,
    height = 1800,
    res = 220
  )

  pheatmap(
    tissue_dist_mat,
    clustering_distance_rows = tissue_dists,
    clustering_distance_cols = tissue_dists,
    annotation_col = tissue_ann,
    annotation_row = tissue_ann,
    annotation_colors = list(dpi = ann_colors$dpi),
    fontsize = 8,
    main = paste0("Sample distance heatmap: ", tissue_name)
  )

  dev.off()

  normalized_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene")

  write.csv(
    normalized_counts,
    file.path(TABLE_DIR, paste0("GSE310471_", tissue_name, "_normalized_counts.csv")),
    row.names = FALSE
  )

  list(
    dds = dds,
    vst = vst_obj,
    vst_mat = vst_mat,
    sample_info = tissue_samples
  )
}

tissue_results <- list(
  Lung = run_tissue_deseq("Lung"),
  Tonsil = run_tissue_deseq("Tonsil")
)

filter_summary <- purrr::map_dfr(names(tissue_results), function(tissue_name) {
  tibble::tibble(
    tissue = tissue_name,
    input_genes = nrow(count_mat),
    min_count = MIN_COUNT,
    min_samples = MIN_SAMPLES,
    genes_used_for_deseq2 = nrow(tissue_results[[tissue_name]]$dds)
  )
})

write.csv(
  filter_summary,
  file.path(TABLE_DIR, "GSE310471_tissue_filter_summary.csv"),
  row.names = FALSE
)

# ------------------------- 7. Contrast results --------------------------------

contrast_plan <- tibble::tribble(
  ~contrast_name, ~tissue,  ~dpi_level,
  "Lung_3DPI_vs_baseline",   "Lung",   "3DPI",
  "Lung_4DPI_vs_baseline",   "Lung",   "4DPI",
  "Lung_5DPI_vs_baseline",   "Lung",   "5DPI",
  "Tonsil_3DPI_vs_baseline", "Tonsil", "3DPI",
  "Tonsil_4DPI_vs_baseline", "Tonsil", "4DPI",
  "Tonsil_5DPI_vs_baseline", "Tonsil", "5DPI"
)

map_human_annotations <- function(genes) {

  symbols <- genes

  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(symbols),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  tibble::tibble(
    Gene = symbols,
    human_ENTREZID = unname(entrez[symbols])
  )
}

plot_volcano <- function(res_df, contrast_name) {

  plot_df <- res_df %>%
    dplyr::mutate(
      padj_plot = ifelse(is.na(padj), 1, padj),
      neg_log10_padj = -log10(padj_plot),
      direction = dplyr::case_when(
        padj < ALPHA & log2FoldChange >= LOG2FC_CUTOFF ~ "Up",
        padj < ALPHA & log2FoldChange <= -LOG2FC_CUTOFF ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  top_labels <- plot_df %>%
    dplyr::filter(direction != "Not significant") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 15)

  p <- ggplot(plot_df, aes(log2FoldChange, neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.65, size = 1.1) +
    geom_vline(
      xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_hline(
      yintercept = -log10(ALPHA),
      linetype = "dashed",
      color = "grey45"
    ) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = Gene),
      size = 3,
      max.overlaps = 30,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Up" = "#D95F59",
        "Down" = "#2A9D8F",
        "Not significant" = "grey75"
      )
    ) +
    labs(
      title = paste0("GSE310471: ", contrast_name),
      x = "log2 fold-change",
      y = "-log10 adjusted p value",
      color = ""
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(FIG_DIR, paste0("Volcano_", contrast_name, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

plot_top_heatmap <- function(vst_mat, res_df, sample_info_subset, contrast_name, n = 50) {

  top_genes <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = min(n, nrow(.))) %>%
    dplyr::pull(Gene)

  if (length(top_genes) < 2) return(invisible(NULL))

  mat <- vst_mat[top_genes, , drop = FALSE]

  ann <- sample_info_subset %>%
    dplyr::select(sample_id, dpi) %>%
    tibble::column_to_rownames("sample_id")

  png(
    file.path(FIG_DIR, paste0("Heatmap_top", n, "_", contrast_name, ".png")),
    width = 2000,
    height = 2600,
    res = 220
  )

  pheatmap(
    mat,
    scale = "row",
    annotation_col = ann,
    annotation_colors = list(dpi = ann_colors$dpi),
    show_colnames = TRUE,
    fontsize_row = 7,
    main = paste0("Top ", n, " genes: ", contrast_name)
  )

  dev.off()
}

run_go <- function(deg_df, background_genes, contrast_name) {

  if (!has_clusterProfiler) return(NULL)

  suppressPackageStartupMessages(library(clusterProfiler))

  deg_entrez <- deg_df %>%
    dplyr::filter(!is.na(human_ENTREZID), human_ENTREZID != "") %>%
    dplyr::pull(human_ENTREZID) %>%
    unique()

  bg_entrez <- map_human_annotations(background_genes) %>%
    dplyr::filter(!is.na(human_ENTREZID), human_ENTREZID != "") %>%
    dplyr::pull(human_ENTREZID) %>%
    unique()

  if (length(deg_entrez) < 10) {
    message(contrast_name, ": fewer than 10 mapped DEG Entrez IDs; skipping GO.")
    return(NULL)
  }

  ego <- enrichGO(
    gene = deg_entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  go_table <- as.data.frame(ego)

  write.csv(
    go_table,
    file.path(TABLE_DIR, paste0("GSE310471_GO_BP_", contrast_name, ".csv")),
    row.names = FALSE
  )

  if (nrow(go_table) > 0) {

    p_go <- dotplot(ego, showCategory = min(20, nrow(go_table))) +
      ggtitle(paste0("GO Biological Process: ", contrast_name)) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold")
      )

    ggsave(
      file.path(FIG_DIR, paste0("GO_BP_dotplot_", contrast_name, ".png")),
      p_go,
      width = 9,
      height = 7,
      dpi = 300
    )
  }

  ego
}

all_contrast_results <- list()
deg_results <- list()
go_results <- list()

for (i in seq_len(nrow(contrast_plan))) {

  contrast_name <- contrast_plan$contrast_name[i]
  tissue_name <- contrast_plan$tissue[i]
  dpi_level <- contrast_plan$dpi_level[i]

  dds <- tissue_results[[tissue_name]]$dds
  vst_mat <- tissue_results[[tissue_name]]$vst_mat
  sample_info_subset <- tissue_results[[tissue_name]]$sample_info

  res <- results(dds, contrast = c("dpi", dpi_level, "baseline"), alpha = ALPHA)

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::arrange(padj)

  anno <- map_human_annotations(res_df$Gene)

  res_df <- res_df %>%
    dplyr::left_join(anno, by = "Gene")

  deg_df <- res_df %>%
    dplyr::filter(
      !is.na(padj),
      padj < ALPHA,
      abs(log2FoldChange) >= LOG2FC_CUTOFF
    )

  all_contrast_results[[contrast_name]] <- res_df
  deg_results[[contrast_name]] <- deg_df

  write.csv(
    res_df,
    file.path(TABLE_DIR, paste0("GSE310471_DESeq2_all_", contrast_name, ".csv")),
    row.names = FALSE
  )

  write.csv(
    deg_df,
    file.path(TABLE_DIR, paste0(
      "GSE310471_DEG_",
      contrast_name,
      "_padj_0.05_log2FC_1.csv"
    )),
    row.names = FALSE
  )

  write.csv(
    res_df %>% dplyr::slice_head(n = 100),
    file.path(TABLE_DIR, paste0("GSE310471_top100_", contrast_name, ".csv")),
    row.names = FALSE
  )

  plot_volcano(res_df, contrast_name)
  plot_top_heatmap(vst_mat, res_df, sample_info_subset, contrast_name, n = 50)

  go_results[[contrast_name]] <- run_go(deg_df, rownames(dds), contrast_name)

  message(
    contrast_name,
    ": ",
    nrow(deg_df),
    " DE genes at padj < ",
    ALPHA,
    " and |log2FC| >= ",
    LOG2FC_CUTOFF
  )
}

deg_summary <- tibble::tibble(
  contrast = names(deg_results),
  n_deg = vapply(deg_results, nrow, integer(1)),
  n_up = vapply(
    deg_results,
    function(x) sum(x$log2FoldChange >= LOG2FC_CUTOFF),
    integer(1)
  ),
  n_down = vapply(
    deg_results,
    function(x) sum(x$log2FoldChange <= -LOG2FC_CUTOFF),
    integer(1)
  )
)

write.csv(
  deg_summary,
  file.path(TABLE_DIR, "GSE310471_DEG_summary.csv"),
  row.names = FALSE
)

# ---------------------- 8. HUVEC signature tracking ---------------------------

huvec_signature <- c(
  "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2", "IFIT3",
  "IFIT5", "IFIH1", "DDX58", "IRF7", "STAT1", "STAT2", "RSAD2", "ISG15",
  "CXCL10", "CXCL11", "CCL5", "HERC5", "HERC6", "SAMD9", "SAMD9L", "PARP9",
  "USP18", "TRIM21", "ICAM1", "VCAM1", "SELE", "IL6", "ANGPT2"
)

signature_summary <- purrr::map_dfr(names(all_contrast_results), function(contrast_name) {

  all_contrast_results[[contrast_name]] %>%
    dplyr::filter(Gene %in% huvec_signature) %>%
    dplyr::mutate(contrast = contrast_name) %>%
    dplyr::select(
      contrast,
      Gene,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )
})

write.csv(
  signature_summary,
  file.path(TABLE_DIR, "GSE310471_HUVEC_signature_contrast_summary.csv"),
  row.names = FALSE
)

plot_signature_heatmap <- function(tissue_name) {

  vst_mat <- tissue_results[[tissue_name]]$vst_mat
  sample_info_subset <- tissue_results[[tissue_name]]$sample_info

  genes_present <- intersect(huvec_signature, rownames(vst_mat))

  if (length(genes_present) < 2) return(invisible(NULL))

  mat <- vst_mat[genes_present, , drop = FALSE]

  ann <- sample_info_subset %>%
    dplyr::select(sample_id, dpi) %>%
    tibble::column_to_rownames("sample_id")

  png(
    file.path(FIG_DIR, paste0("Heatmap_HUVEC_signature_", tissue_name, ".png")),
    width = 2400,
    height = 3000,
    res = 300
  )

  pheatmap(
    mat,
    scale = "row",
    annotation_col = ann,
    annotation_colors = list(dpi = ann_colors$dpi),
    fontsize = 11,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    cellheight = 18,
    cellwidth = 24
  )

  dev.off()
}

plot_signature_heatmap("Lung")
plot_signature_heatmap("Tonsil")

# ----------------------------- 9. Save objects --------------------------------

saveRDS(
  list(
    sample_info = sample_info,
    count_mat = count_mat,
    tissue_results = tissue_results,
    contrast_plan = contrast_plan,
    all_contrast_results = all_contrast_results,
    deg_results = deg_results,
    deg_summary = deg_summary,
    go_results = go_results,
    huvec_signature = huvec_signature,
    signature_summary = signature_summary
  ),
  file.path(OUT_DIR, "GSE310471_analysis_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

message("Analysis complete.")
message("Results folder: ", OUT_DIR)
