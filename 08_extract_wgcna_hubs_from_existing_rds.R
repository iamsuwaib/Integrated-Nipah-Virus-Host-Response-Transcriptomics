###############################################################################
# Extract WGCNA hub genes from existing WGCNA RDS objects
#
# Use this after 03_wgcna_gse310471_exploratory.R has already completed.
# It avoids rerunning WGCNA and only rebuilds the hub-gene tables.
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "WGCNA")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "WGCNA") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

library(tidyverse)
library(WGCNA)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
table_dir <- file.path(project_dir, "advanced_analyses", "tables")

extract_hubs <- function(tissue_name) {
  rds_path <- file.path(table_dir, paste0("WGCNA_", tissue_name, "_object.rds"))
  if (!file.exists(rds_path)) stop("Missing RDS object: ", rds_path)

  obj <- readRDS(rds_path)
  datExpr <- obj$datExpr
  net <- obj$net
  MEs <- obj$MEs

  gene_module_table <- tibble(
    gene = colnames(datExpr),
    module_color = labels2colors(net$colors),
    module_numeric = net$colors,
    tissue = tissue_name
  )

  hub_table <- map_dfr(sort(unique(net$colors)), function(mod_num) {
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

  readr::write_csv(
    hub_table,
    file.path(table_dir, paste0("WGCNA_", tissue_name, "_hub_genes_top25_per_module.csv"))
  )

  hub_table
}

lung_hubs <- extract_hubs("Lung")
tonsil_hubs <- extract_hubs("Tonsil")

combined <- bind_rows(lung_hubs, tonsil_hubs)
readr::write_csv(
  combined,
  file.path(table_dir, "WGCNA_combined_hub_genes_top25_per_module.csv")
)

message("Hub extraction complete.")
message("Rows written: ", nrow(combined))

