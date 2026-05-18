###############################################################################
# Exploratory WGCNA / module analysis for GSE310471
#
# Important:
#   WGCNA is most powerful with larger sample sizes. Here, lung and tonsil have
#   modest n, so module results should be treated as exploratory and used to
#   support the integrated signature story, not as standalone proof.
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "DESeq2", "WGCNA", "pheatmap", "RColorBrewer")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("DESeq2", "WGCNA")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

library(tidyverse)
library(DESeq2)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)

allowWGCNAThreads()

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(project_dir, "GSE310471/GSE310471_counts.tsv/GSE310471_counts.tsv")
metadata_file <- file.path(project_dir, "GSE310471/results/tables/GSE310471_sample_metadata.csv")

counts <- readr::read_tsv(count_file, show_col_types = FALSE)
metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)

# Some processed count matrices contain duplicated gene symbols. WGCNA requires
# unique row names, so duplicate symbols are collapsed by summing counts before
# the expression matrix is created.
duplicate_gene_summary <- counts %>%
  count(Gene, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  arrange(desc(n_rows), Gene)

readr::write_csv(
  duplicate_gene_summary,
  file.path(table_dir, "GSE310471_duplicate_gene_rows_collapsed_for_WGCNA.csv")
)

counts_collapsed <- counts %>%
  group_by(Gene) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

count_mat <- counts_collapsed %>%
  column_to_rownames("Gene") %>%
  as.matrix()
storage.mode(count_mat) <- "integer"

metadata <- metadata %>%
  filter(sample_id %in% colnames(count_mat)) %>%
  mutate(
    dpi_numeric = case_when(
      dpi == "baseline" ~ 0,
      TRUE ~ as.numeric(str_extract(dpi, "[0-9]+"))
    ),
    infected = ifelse(dpi == "baseline", 0, 1),
    dpi_factor = factor(dpi, levels = c("baseline", "3DPI", "4DPI", "5DPI"))
  ) %>%
  arrange(match(sample_id, colnames(count_mat)))

run_wgcna_for_tissue <- function(tissue_name, top_variable_genes = 8000) {
  message("Running WGCNA for ", tissue_name)

  meta_t <- metadata %>% filter(tissue == tissue_name)
  mat_t <- count_mat[, meta_t$sample_id, drop = FALSE]

  keep <- rowSums(mat_t >= 10) >= max(3, floor(ncol(mat_t) * 0.25))
  mat_t <- mat_t[keep, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = mat_t,
    colData = as.data.frame(meta_t) %>% column_to_rownames("sample_id"),
    design = ~ dpi_factor
  )
  dds <- estimateSizeFactors(dds)
  vst_mat <- assay(vst(dds, blind = TRUE))

  gene_var <- apply(vst_mat, 1, var, na.rm = TRUE)
  selected_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(top_variable_genes, length(gene_var)))]
  datExpr <- t(vst_mat[selected_genes, , drop = FALSE])

  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }

  powers <- c(1:10, seq(12, 20, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
  soft_power <- sft$powerEstimate
  if (is.na(soft_power)) soft_power <- 6

  png(file.path(figure_dir, paste0("WGCNA_", tissue_name, "_soft_threshold.png")),
      width = 2400, height = 1200, res = 300)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft threshold power", ylab = "Scale-free topology model fit",
       type = "n", main = paste(tissue_name, "scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, col = "firebrick")
  abline(h = 0.8, col = "grey50", lty = 2)
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft threshold power", ylab = "Mean connectivity",
       type = "n", main = paste(tissue_name, "mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "firebrick")
  dev.off()

  net <- blockwiseModules(
    datExpr,
    power = soft_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3
  )

  module_colors <- labels2colors(net$colors)
  MEs <- orderMEs(net$MEs)

  trait_df <- meta_t %>%
    filter(sample_id %in% rownames(datExpr)) %>%
    arrange(match(sample_id, rownames(datExpr))) %>%
    transmute(
      dpi_numeric = dpi_numeric,
      infected = infected,
      baseline = as.integer(dpi == "baseline"),
      dpi3 = as.integer(dpi == "3DPI"),
      dpi4 = as.integer(dpi == "4DPI"),
      dpi5 = as.integer(dpi == "5DPI")
    )

  module_trait_cor <- cor(MEs, trait_df, use = "p")
  module_trait_p <- corPvalueStudent(module_trait_cor, nrow(datExpr))

  cor_table <- as.data.frame(module_trait_cor) %>%
    rownames_to_column("module") %>%
    pivot_longer(-module, names_to = "trait", values_to = "correlation") %>%
    left_join(
      as.data.frame(module_trait_p) %>%
        rownames_to_column("module") %>%
        pivot_longer(-module, names_to = "trait", values_to = "pvalue"),
      by = c("module", "trait")
    ) %>%
    mutate(tissue = tissue_name, soft_power = soft_power) %>%
    arrange(pvalue)

  gene_module_table <- tibble(
    gene = colnames(datExpr),
    module_color = module_colors,
    module_numeric = net$colors
  ) %>%
    mutate(tissue = tissue_name)

  hub_tables <- map_dfr(sort(unique(net$colors)), function(mod_num) {
    mod_info <- gene_module_table %>% filter(module_numeric == mod_num)
    mod_genes <- mod_info$gene
    if (length(mod_genes) < 5) return(NULL)
    ME_name <- paste0("ME", mod_num)
    if (!ME_name %in% colnames(MEs)) return(NULL)
    kME <- cor(datExpr[, mod_genes, drop = FALSE], MEs[, ME_name], use = "p")
    tibble(
      tissue = tissue_name,
      module_color = unique(mod_info$module_color)[1],
      module_numeric = mod_num,
      gene = mod_genes,
      kME = as.numeric(kME[, 1])
    ) %>%
      arrange(desc(abs(kME))) %>%
      slice_head(n = 25)
  })

  readr::write_csv(cor_table, file.path(table_dir, paste0("WGCNA_", tissue_name, "_module_trait_correlations.csv")))
  readr::write_csv(gene_module_table, file.path(table_dir, paste0("WGCNA_", tissue_name, "_gene_modules.csv")))
  readr::write_csv(hub_tables, file.path(table_dir, paste0("WGCNA_", tissue_name, "_hub_genes_top25_per_module.csv")))

  png(file.path(figure_dir, paste0("WGCNA_", tissue_name, "_module_trait_heatmap.png")),
      width = 2200, height = 1800, res = 300)
  pheatmap(
    module_trait_cor,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
    breaks = seq(-1, 1, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = matrix(sprintf("%.2f\np=%.2g", module_trait_cor, module_trait_p),
                             nrow = nrow(module_trait_cor)),
    fontsize = 8,
    main = paste0("WGCNA module-trait correlations: ", tissue_name)
  )
  dev.off()

  saveRDS(
    list(
      tissue = tissue_name,
      datExpr = datExpr,
      net = net,
      MEs = MEs,
      traits = trait_df,
      module_trait_cor = module_trait_cor,
      module_trait_p = module_trait_p,
      soft_power = soft_power
    ),
    file.path(table_dir, paste0("WGCNA_", tissue_name, "_object.rds"))
  )

  invisible(TRUE)
}

run_wgcna_for_tissue("Lung", top_variable_genes = 8000)
run_wgcna_for_tissue("Tonsil", top_variable_genes = 8000)

message("Exploratory WGCNA complete.")
