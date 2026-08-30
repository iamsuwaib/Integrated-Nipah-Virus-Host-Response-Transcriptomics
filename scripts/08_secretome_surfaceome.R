# ============================================================
# 08_secretome_surfaceome.R
# Proteomics-relevant extracellular and cell-surface mediator prioritization for NiV transcriptomics
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(forcats)
})

# -----------------------------
# 0. Paths
# -----------------------------

ROOT <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"

integrated_file <- file.path(ROOT, "integrated_results/tables/integrated_signature_long.csv")
biomarker_file  <- file.path(ROOT, "advanced_analyses/tables/candidate_biomarker_table_full.csv")

out_table_dir <- file.path(ROOT, "advanced_analyses/tables")
out_fig_dir   <- file.path(ROOT, "advanced_analyses/figures")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Helper functions
# -----------------------------

check_required_cols <- function(df, required_cols, file_name) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", file_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max_abs <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  max(abs(x))
}

# -----------------------------
# 2. Load data
# -----------------------------

integrated <- read_csv(integrated_file, show_col_types = FALSE)

check_required_cols(
  integrated,
  c("gene", "dataset", "tissue", "timepoint", "log2FC", "padj", "significant"),
  "integrated_signature_long.csv"
)

integrated <- integrated %>%
  mutate(
    gene = toupper(gene),
    significant = as.logical(significant),
    direction = case_when(
      is.na(log2FC) ~ NA_character_,
      log2FC > 0 ~ "up",
      log2FC < 0 ~ "down",
      TRUE ~ "neutral"
    )
  )

if (file.exists(biomarker_file)) {
  biomarkers <- read_csv(biomarker_file, show_col_types = FALSE) %>%
    mutate(gene = toupper(gene))
} else {
  warning("Biomarker file not found. Continuing without biomarker annotation.")
  biomarkers <- tibble(gene = character())
}

# -----------------------------
# 3. Curated secretome/surfaceome knowledge
# -----------------------------

chemokines <- c(
  paste0("CCL", 1:28),
  paste0("CXCL", 1:17),
  "CX3CL1",
  paste0("XCL", 1:2)
)

cytokines <- c(
  paste0("IL", 1:38),
  "IL1A", "IL1B", "IL1RN",
  "IL6", "IL10", "IL12A", "IL12B", "IL15", "IL18", "IL33",
  "TNF", "TNFSF10", "TNFSF13B", "TNFSF14",
  "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8",
  "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21",
  "IFNB1", "IFNG", "IFNL1", "IFNL2", "IFNL3",
  "CSF1", "CSF2", "CSF3",
  "TGFB1", "TGFB2", "TGFB3"
)

complement <- c(
  "C1QA", "C1QB", "C1QC", "C1R", "C1S",
  "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  "CFB", "CFD", "CFH", "CFI", "CFP",
  "CD46", "CD55", "CD59",
  "C3AR1", "C5AR1", "C5AR2"
)

coagulation_fibrinolysis <- c(
  "F2", "F3", "F5", "F7", "F8", "F9", "F10", "F11", "F12",
  "F13A1", "F13B",
  "VWF", "SERPINE1", "SERPINC1", "SERPINA1", "THBD",
  "PLAT", "PLAU", "PLAUR", "PROS1", "PROC",
  "FGA", "FGB", "FGG"
)

ecm_remodeling <- c(
  paste0("MMP", 1:28),
  paste0("TIMP", 1:4),
  "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
  "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL8A2",
  "FN1", "SPP1", "THBS1", "THBS2", "VCAN",
  "LAMB1", "LAMB2", "LAMB3",
  "LAMC1", "LAMC2",
  "ELN", "FBLN1", "FBLN2",
  "POSTN", "LOX", "LOXL1", "LOXL2",
  "HAS1", "HAS2", "HAS3"
)

adhesion_endothelial <- c(
  "ICAM1", "ICAM2", "VCAM1",
  "SELE", "SELL", "SELP", "SELPLG",
  "PECAM1", "MCAM", "MADCAM1", "ALCAM", "NCAM1",
  "CDH5", "CDH1", "CDH2",
  "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6",
  "ITGAL", "ITGAM", "ITGAX",
  "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8"
)

receptors_surface <- c(
  paste0("TLR", 1:10),
  "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2",
  "IL1R1", "IL1R2", "IL2RA", "IL2RB", "IL2RG", "IL4R",
  "IL6R", "IL6ST", "IL7R", "IL10RA", "IL10RB",
  "IL12RB1", "IL12RB2", "IL13RA1", "IL15RA",
  "IL17RA", "IL18R1", "IL18RAP",
  paste0("CCR", 1:10),
  paste0("CXCR", 1:6),
  "CX3CR1",
  "TNFRSF1A", "TNFRSF1B", "TNFRSF4", "TNFRSF8", "TNFRSF9",
  "TNFRSF10A", "TNFRSF10B", "TNFRSF11A", "TNFRSF13B",
  "TNFRSF14", "TNFRSF17",
  "KDR", "FLT1", "FLT4", "TEK", "ENG",
  "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
  "TREM1", "TREM2",
  "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B",
  "SIGLEC1", "SIGLEC5", "SIGLEC7", "SIGLEC9"
)

antigen_presentation_surface <- c(
  "HLA-A", "HLA-B", "HLA-C",
  "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
  "HLA-DPA1", "HLA-DPB1",
  "B2M", "CD74", "TAP1", "TAP2"
)

acute_phase_secreted <- c(
  "S100A8", "S100A9", "S100A12",
  "LCN2", "HP", "HPX", "APOE", "APOA1", "APOB",
  "VEGFA", "VEGFB", "PGF", "ANGPT1", "ANGPT2"
)

cd_surface <- paste0("CD", 1:400)

secretome_genes <- unique(c(
  chemokines,
  cytokines,
  complement,
  coagulation_fibrinolysis,
  ecm_remodeling,
  acute_phase_secreted
))

surfaceome_genes <- unique(c(
  adhesion_endothelial,
  receptors_surface,
  antigen_presentation_surface,
  cd_surface,
  complement
))

# -----------------------------
# 4. Annotation functions
# -----------------------------

annotate_compartment <- function(gene) {
  secreted <- gene %in% secretome_genes
  
  surface <- gene %in% surfaceome_genes ||
    str_detect(gene, "^CD[0-9]+$") ||
    str_detect(gene, "^HLA-") ||
    str_detect(gene, "^ITGA|^ITGB|^TLR|^CCR|^CXCR|^TNFRSF")
  
  case_when(
    secreted & surface ~ "secreted_or_surface-associated",
    secreted ~ "secretome",
    surface ~ "surfaceome",
    TRUE ~ "other"
  )
}

annotate_physiology <- function(gene) {
  case_when(
    gene %in% chemokines ~ "chemokine signaling",
    gene %in% cytokines ~ "cytokine/interferon mediator",
    gene %in% complement ~ "complement cascade",
    gene %in% coagulation_fibrinolysis ~ "coagulation/fibrinolysis",
    gene %in% ecm_remodeling ~ "extracellular matrix/remodeling",
    gene %in% adhesion_endothelial ~ "vascular adhesion/endothelial activation",
    gene %in% antigen_presentation_surface ~ "antigen presentation",
    gene %in% receptors_surface ~ "surface receptor/signaling",
    str_detect(gene, "^CD[0-9]+$") ~ "cell-surface immune marker",
    gene %in% acute_phase_secreted ~ "acute-phase/vascular secreted mediator",
    TRUE ~ "other extracellular/surface candidate"
  )
}

# -----------------------------
# 5. Candidate table
# -----------------------------

candidate_long <- integrated %>%
  mutate(
    compartment = map_chr(gene, annotate_compartment),
    physiology_class = map_chr(gene, annotate_physiology),
    context_id = paste(dataset, tissue, timepoint, sep = " | ")
  ) %>%
  filter(compartment != "other")

candidate_summary <- candidate_long %>%
  filter(scoring_included) %>%
  group_by(gene, compartment, physiology_class) %>%
  summarise(
    n_detected = n(),
    n_significant = sum(significant, na.rm = TRUE),
    n_significant_up = sum(significant & direction == "up", na.rm = TRUE),
    n_significant_down = sum(significant & direction == "down", na.rm = TRUE),
    n_datasets = n_distinct(dataset),
    n_tissues = n_distinct(tissue),
    n_timepoints = n_distinct(timepoint),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    median_log2FC = median(log2FC, na.rm = TRUE),
    max_abs_log2FC = safe_max_abs(log2FC),
    min_padj = safe_min(padj),
    detected_contexts = paste(unique(context_id), collapse = "; "),
    significant_up_contexts = paste(unique(context_id[significant & direction == "up"]), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    direction_bias = case_when(
      n_significant_up > n_significant_down ~ "mostly_up",
      n_significant_down > n_significant_up ~ "mostly_down",
      n_significant_up == n_significant_down & n_significant_up > 0 ~ "mixed",
      TRUE ~ "not_significant"
    ),
    priority_score =
      (2.5 * n_significant_up) +
      (1.5 * n_datasets) +
      n_tissues +
      n_timepoints +
      pmin(replace_na(max_abs_log2FC, 0), 5) -
      (1.0 * n_significant_down),
    evidence_level = case_when(
      n_significant_up >= 6 & n_datasets >= 2 ~ "high",
      n_significant_up >= 3 & n_datasets >= 2 ~ "moderate",
      n_significant_up >= 2 ~ "supportive",
      n_significant_up >= 1 ~ "limited",
      TRUE ~ "detected_not_significant"
    ),
    translational_relevance = case_when(
      compartment == "secretome" ~
        "Secreted mediator; suitable for plasma/serum/secretome validation",
      compartment == "surfaceome" ~
        "Cell-surface mediator; suitable for flow cytometry, CyTOF, spatial or membrane proteomics validation",
      compartment == "secreted_or_surface-associated" ~
        "Extracellular or membrane-associated mediator; suitable for both soluble and cell-surface follow-up",
      TRUE ~ "Candidate extracellular/surface mediator"
    )
  ) %>%
  arrange(desc(priority_score), desc(n_significant_up), min_padj)

# -----------------------------
# 6. Add biomarker annotations, if present
# -----------------------------

optional_biomarker_cols <- intersect(
  c("gene", "manuscript_module", "evidence_strength", "conserved_call", "concise_rationale"),
  colnames(biomarkers)
)

if (length(optional_biomarker_cols) > 1) {
  candidate_summary <- candidate_summary %>%
    left_join(
      biomarkers %>%
        select(all_of(optional_biomarker_cols)) %>%
        distinct(gene, .keep_all = TRUE),
      by = "gene"
    )
}

write_csv(
  candidate_long,
  file.path(out_table_dir, "secretome_surfaceome_candidate_long.csv")
)

write_csv(
  candidate_summary,
  file.path(out_table_dir, "secretome_surfaceome_candidate_table_refined.csv")
)

# -----------------------------
# 7. High-confidence shortlist
# -----------------------------

high_confidence_candidates <- candidate_summary %>%
  filter(
    evidence_level %in% c("high", "moderate", "supportive"),
    n_significant_up >= 2,
    mean_log2FC > 0
  ) %>%
  arrange(desc(priority_score))

write_csv(
  high_confidence_candidates,
  file.path(out_table_dir, "secretome_surfaceome_high_confidence_shortlist.csv")
)

# -----------------------------
# 8. Category summary
# -----------------------------

category_summary <- candidate_summary %>%
  group_by(compartment, physiology_class) %>%
  summarise(
    n_genes = n(),
    n_high = sum(evidence_level == "high"),
    n_moderate = sum(evidence_level == "moderate"),
    n_supportive = sum(evidence_level == "supportive"),
    high_or_moderate = sum(evidence_level %in% c("high", "moderate")),
    median_significant_up = median(n_significant_up, na.rm = TRUE),
    top_genes = paste(head(gene[order(-priority_score)], 10), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(high_or_moderate), desc(n_genes))

write_csv(
  category_summary,
  file.path(out_table_dir, "secretome_surfaceome_category_summary_refined.csv")
)

# -----------------------------
# 9. Top-candidate lollipop plot
# -----------------------------

top_candidates <- candidate_summary %>%
  filter(n_significant_up > 0) %>%
  slice_max(priority_score, n = 30) %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),
    evidence_level = factor(
      evidence_level,
      levels = c("high", "moderate", "supportive", "limited", "detected_not_significant")
    )
  )

p1 <- ggplot(top_candidates, aes(x = priority_score, y = gene, color = compartment)) +
  geom_segment(aes(x = 0, xend = priority_score, yend = gene), linewidth = 0.7, alpha = 0.6) +
  geom_point(aes(size = n_significant_up, shape = evidence_level), alpha = 0.95) +
  scale_color_manual(values = c(
    "secretome" = "#D55E00",
    "surfaceome" = "#0072B2",
    "secreted_or_surface-associated" = "#009E73"
  )) +
  labs(
    title = "Extracellular and cell-surface mediators in NiV host-response signatures",
    subtitle = "Ranking integrates recurrent upregulation, dataset/tissue breadth and effect size",
    x = "Priority score",
    y = "Gene",
    color = "Compartment",
    size = "Significant upregulated contexts",
    shape = "Evidence level"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(out_fig_dir, "secretome_surfaceome_top_candidates_refined.png"),
  p1,
  width = 10,
  height = 8,
  dpi = 300
)

# -----------------------------
# 10. Physiology-class summary plot
# -----------------------------

p2 <- category_summary %>%
  mutate(
    physiology_class = fct_reorder(physiology_class, n_genes)
  ) %>%
  ggplot(aes(x = n_genes, y = physiology_class, fill = compartment)) +
  geom_col(width = 0.75, color = "grey25") +
  labs(
    title = "Extracellular and cell-surface mediator classes",
    x = "Number of candidate genes",
    y = "Physiology class",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "secretome" = "#D55E00",
    "surfaceome" = "#0072B2",
    "secreted_or_surface-associated" = "#009E73"
  )) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(out_fig_dir, "secretome_surfaceome_category_summary_refined.png"),
  p2,
  width = 11,
  height = 6.5,
  dpi = 300
)

# -----------------------------
# 11. Heatmap-style tile plot
# -----------------------------

top_heatmap_genes <- candidate_summary %>%
  filter(n_significant_up > 0) %>%
  slice_max(priority_score, n = 40) %>%
  pull(gene)

heatmap_df <- candidate_long %>%
  filter(gene %in% top_heatmap_genes) %>%
  mutate(
    contrast_label = paste(dataset, tissue, timepoint, sep = " | "),
    gene = factor(gene, levels = rev(top_heatmap_genes)),
    log2FC_capped = pmax(pmin(log2FC, 5), -5),
    significance_label = ifelse(significant, "significant", "not significant")
  )

p3 <- ggplot(heatmap_df, aes(x = contrast_label, y = gene, fill = log2FC_capped)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(
    data = heatmap_df %>% filter(significant),
    aes(x = contrast_label, y = gene),
    size = 0.6,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-5, 5),
    name = "log2FC"
  ) +
  labs(
    title = "Expression of extracellular and cell-surface mediator candidates across NiV transcriptomic contexts",
    subtitle = "Black dots indicate significant contrasts",
    x = "Dataset | tissue/model | time point",
    y = "Candidate gene"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_fig_dir, "secretome_surfaceome_heatmap_refined.png"),
  p3,
  width = 13,
  height = 9,
  dpi = 300
)

# -----------------------------
# 12. Separate secretome and surfaceome shortlists
# -----------------------------

secretome_shortlist <- candidate_summary %>%
  filter(compartment %in% c("secretome", "secreted_or_surface-associated")) %>%
  filter(n_significant_up > 0) %>%
  arrange(desc(priority_score))

surfaceome_shortlist <- candidate_summary %>%
  filter(compartment %in% c("surfaceome", "secreted_or_surface-associated")) %>%
  filter(n_significant_up > 0) %>%
  arrange(desc(priority_score))

write_csv(
  secretome_shortlist,
  file.path(out_table_dir, "secretome_candidate_shortlist_refined.csv")
)

write_csv(
  surfaceome_shortlist,
  file.path(out_table_dir, "surfaceome_candidate_shortlist_refined.csv")
)

# -----------------------------
# 13. Interpretation notes
# -----------------------------

top20 <- candidate_summary %>%
  slice_head(n = 20) %>%
  pull(gene)

top_secreted <- secretome_shortlist %>%
  slice_head(n = 15) %>%
  pull(gene)

top_surface <- surfaceome_shortlist %>%
  slice_head(n = 15) %>%
  pull(gene)

interpretation <- c(
  "Proteomics-relevant extracellular/cell-surface mediator analysis summary",
  "====================================",
  paste0("Total extracellular/surface candidates detected: ", nrow(candidate_summary)),
  paste0(
    "High/moderate/supportive candidates: ",
    sum(candidate_summary$evidence_level %in% c("high", "moderate", "supportive"))
  ),
  "",
  "Top overall candidates:",
  paste(top20, collapse = ", "),
  "",
  "Top secretome candidates:",
  paste(top_secreted, collapse = ", "),
  "",
  "Top surfaceome candidates:",
  paste(top_surface, collapse = ", "),
  "",
  "Main physiology classes:",
  paste(
    category_summary$physiology_class,
    category_summary$n_genes,
    sep = " = ",
    collapse = "; "
  ),
  "",
  "Manuscript interpretation:",
  "This analysis links conserved NiV transcriptional programs to proteomics-relevant extracellular and cell-surface mediators of tissue-level disease physiology.",
  "Extracellular mediator candidates highlight soluble inflammatory, chemokine, complement, coagulation, acute-phase, and matrix-remodeling pathways.",
  "Cell-surface mediator candidates highlight endothelial activation, immune-cell interaction, antigen-presentation, and receptor-signaling components.",
  "These candidates provide a prioritized shortlist for plasma proteomics, secretome profiling, membrane proteomics, flow cytometry, spatial validation, or targeted mechanistic follow-up."
)

writeLines(
  interpretation,
  file.path(out_table_dir, "secretome_surfaceome_interpretation_notes_refined.txt")
)

# -----------------------------
# 14. Final messages
# -----------------------------

message("Secretome/surfaceome refined analysis complete.")
message("Candidate table: ", file.path(out_table_dir, "secretome_surfaceome_candidate_table_refined.csv"))
message("High-confidence shortlist: ", file.path(out_table_dir, "secretome_surfaceome_high_confidence_shortlist.csv"))
message("Secretome shortlist: ", file.path(out_table_dir, "secretome_candidate_shortlist_refined.csv"))
message("Surfaceome shortlist: ", file.path(out_table_dir, "surfaceome_candidate_shortlist_refined.csv"))
message("Figures written to: ", out_fig_dir)


