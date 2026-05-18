###############################################################################
# Candidate biomarker table
#
# Goal:
#   Build a clean manuscript-ready candidate table separating:
#     1. conserved antiviral core
#     2. mainly in vivo complement/coagulation module
#     3. mainly HUVEC/endothelial module
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "openxlsx", "ggplot2", "scales")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}

library(tidyverse)
library(openxlsx)
library(ggplot2)
library(scales)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
integrated_dir <- file.path(project_dir, "integrated_results")
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

top_conserved_path <- file.path(integrated_dir, "tables", "integrated_signature_top_conserved_genes.csv")
long_path <- file.path(integrated_dir, "tables", "integrated_signature_long.csv")

top_conserved <- readr::read_csv(top_conserved_path, show_col_types = FALSE)
signature_long <- readr::read_csv(long_path, show_col_types = FALSE)

candidate_table <- top_conserved %>%
  mutate(
    evidence_strength = case_when(
      conserved_call == "conserved_HUVEC_and_in_vivo" & n_significant_up >= 7 ~ "high",
      conserved_call == "conserved_HUVEC_and_in_vivo" & n_significant_up >= 5 ~ "moderate_high",
      conserved_call == "mainly_in_vivo" & n_in_vivo_significant_up >= 3 ~ "in_vivo_specific",
      conserved_call == "mainly_HUVEC" & n_huvec_significant_up >= 2 ~ "HUVEC_enriched",
      TRUE ~ "supporting"
    ),
    manuscript_module = case_when(
      conserved_call == "conserved_HUVEC_and_in_vivo" &
        category %in% c("ISG_antiviral", "cytokine_chemokine") ~
        "Conserved antiviral/IFN core",
      conserved_call == "mainly_in_vivo" &
        category %in% c("complement", "coagulation") ~
        "In vivo complement/coagulation disease module",
      conserved_call == "mainly_in_vivo" ~
        "Mainly in vivo immune/progression module",
      conserved_call == "mainly_HUVEC" ~
        "HUVEC-enriched endothelial/early-response module",
      TRUE ~ "Supporting candidate"
    ),
    priority_score = (n_significant_up * 2) +
      n_huvec_significant_up +
      n_in_vivo_significant_up +
      ifelse(category %in% c("ISG_antiviral", "cytokine_chemokine"), 2, 0) +
      ifelse(category %in% c("complement", "coagulation"), 1.5, 0),
    concise_rationale = case_when(
      manuscript_module == "Conserved antiviral/IFN core" ~
        "Repeatedly upregulated across HUVEC and in vivo contrasts; supports conserved Nipah antiviral host-response signature.",
      manuscript_module == "In vivo complement/coagulation disease module" ~
        "Predominantly induced in AGM tissues; supports disease-progression biology not captured by HUVEC alone.",
      manuscript_module == "Mainly in vivo immune/progression module" ~
        "Stronger in AGM tissue contrasts; may reflect immune-cell or tissue-level infection progression.",
      manuscript_module == "HUVEC-enriched endothelial/early-response module" ~
        "Supported primarily by human endothelial-cell contrasts; useful for early endothelial response framing.",
      TRUE ~ "Secondary supporting gene for pathway-level interpretation."
    )
  ) %>%
  arrange(desc(priority_score), min_padj) %>%
  select(
    gene, category, manuscript_module, evidence_strength,
    n_detected, n_significant_up, n_huvec_significant_up,
    n_in_vivo_significant_up, mean_log2FC, max_abs_log2FC,
    min_padj, priority_score, conserved_call, concise_rationale
  )

shortlist <- candidate_table %>%
  filter(
    manuscript_module %in% c(
      "Conserved antiviral/IFN core",
      "In vivo complement/coagulation disease module",
      "HUVEC-enriched endothelial/early-response module"
    ),
    evidence_strength %in% c("high", "moderate_high", "in_vivo_specific", "HUVEC_enriched")
  ) %>%
  group_by(manuscript_module) %>%
  slice_max(priority_score, n = 15, with_ties = FALSE) %>%
  ungroup()

readr::write_csv(candidate_table, file.path(table_dir, "candidate_biomarker_table_full.csv"))
readr::write_csv(shortlist, file.path(table_dir, "candidate_biomarker_shortlist.csv"))

wb <- createWorkbook()
addWorksheet(wb, "full_candidate_table")
writeData(wb, "full_candidate_table", candidate_table)
addWorksheet(wb, "shortlist")
writeData(wb, "shortlist", shortlist)
freezePane(wb, "full_candidate_table", firstRow = TRUE)
freezePane(wb, "shortlist", firstRow = TRUE)
saveWorkbook(wb, file.path(table_dir, "candidate_biomarker_table.xlsx"), overwrite = TRUE)

plot_df <- shortlist %>%
  mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    manuscript_module = factor(
      manuscript_module,
      levels = c(
        "Conserved antiviral/IFN core",
        "In vivo complement/coagulation disease module",
        "HUVEC-enriched endothelial/early-response module",
        "Mainly in vivo immune/progression module",
        "Supporting candidate"
      )
    )
  )

png(
  file.path(figure_dir, "candidate_biomarker_shortlist_lollipop.png"),
  width = 3200, height = 2400, res = 300
)
print(
  ggplot(plot_df, aes(x = priority_score, y = gene)) +
    geom_segment(aes(x = 0, xend = priority_score, yend = gene),
                 color = "grey70", linewidth = 0.6) +
    geom_point(
      aes(fill = manuscript_module, size = n_in_vivo_significant_up),
      shape = 21, color = "grey15", stroke = 0.25, alpha = 0.95
    ) +
    scale_fill_manual(
      values = c(
        "Conserved antiviral/IFN core" = "#762A83",
        "In vivo complement/coagulation disease module" = "#E66101",
        "HUVEC-enriched endothelial/early-response module" = "#1F78B4",
        "Mainly in vivo immune/progression module" = "#5AAE61",
        "Supporting candidate" = "grey65"
      )
    ) +
    scale_size_continuous(range = c(2, 8), breaks = c(0, 2, 4, 6)) +
    labs(
      title = "Candidate Nipah host-response biomarkers",
      subtitle = "Ranked using cross-dataset significance, model conservation, and biological category",
      x = "Priority score",
      y = NULL,
      fill = "Module",
      size = "In vivo\nsupport"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "bottom"
    )
)
dev.off()

message("Candidate biomarker outputs written to: ", table_dir)

