###############################################################################
# Focused WGCNA discovery figure for Nipah host-response manuscript
#
# Goal:
#   Convert the candidate WGCNA membership result into a clean manuscript-level
#   figure and table focused only on the two discoveries supported by the data:
#
#   1. Lung antiviral/IFN candidates are strong members of the lung tan module.
#   2. Tonsil complement/coagulation candidates are strong members of the tonsil
#      blue infection/progression module.
#
# Input:
#   advanced_analyses/tables/WGCNA_candidate_signature_membership_summary.csv
#
# Output:
#   advanced_analyses/tables/focused_wgcna_candidate_modules.csv
#   advanced_analyses/tables/focused_wgcna_candidate_modules.xlsx
#   advanced_analyses/figures/focused_wgcna_candidate_modules_panel.png
#   advanced_analyses/figures/focused_wgcna_candidate_modules_heatmap.png
###############################################################################

options(stringsAsFactors = FALSE)

packages <- c("tidyverse", "ggplot2", "openxlsx", "patchwork", "scales")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}

library(tidyverse)
library(ggplot2)
library(openxlsx)
library(patchwork)
library(scales)

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
advanced_dir <- file.path(project_dir, "advanced_analyses")
table_dir <- file.path(advanced_dir, "tables")
figure_dir <- file.path(advanced_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(table_dir, "WGCNA_candidate_signature_membership_summary.csv")
if (!file.exists(input_file)) {
  stop("Missing input file. Run 03c_wgcna_candidate_signature_membership.R first: ", input_file)
}

membership <- readr::read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    gene = as.character(gene),
    tissue = as.character(tissue),
    module_numeric = as.integer(module_numeric),
    kME = as.numeric(kME),
    infected_module_cor = as.numeric(infected_module_cor),
    dpi_module_cor = as.numeric(dpi_module_cor),
    priority_score = as.numeric(priority_score)
  ) %>%
  distinct(tissue, gene, module_numeric, module_eigengene, kME, .keep_all = TRUE)

###############################################################################
# Define the focused biological panels
###############################################################################

lung_antiviral <- membership %>%
  filter(
    tissue == "Lung",
    module_color == "tan",
    module_numeric == 12,
    manuscript_module %in% c(
      "Conserved antiviral/IFN core",
      "Mainly in vivo immune/progression module"
    ),
    strong_module_member,
    strong_infection_module,
    strong_progression_module
  ) %>%
  arrange(desc(priority_score), desc(abs(kME))) %>%
  slice_head(n = 20) %>%
  mutate(
    discovery_panel = "Lung antiviral IFN module",
    module_label = "Lung tan module (ME12)",
    discovery_interpretation =
      "Conserved antiviral/IFN genes form a strong lung infection- and DPI-associated co-expression module."
  )

tonsil_complement <- membership %>%
  filter(
    tissue == "Tonsil",
    module_color == "blue",
    module_numeric == 2,
    manuscript_module == "In vivo complement/coagulation disease module",
    strong_module_member,
    strong_infection_module,
    strong_progression_module
  ) %>%
  arrange(desc(priority_score), desc(abs(kME))) %>%
  mutate(
    discovery_panel = "Tonsil complement/coagulation module",
    module_label = "Tonsil blue module (ME2)",
    discovery_interpretation =
      "Complement and coagulation-associated genes form a strong tonsil infection- and DPI-associated disease module."
  )

focused <- bind_rows(lung_antiviral, tonsil_complement) %>%
  mutate(
    abs_kME = abs(kME),
    panel_order = case_when(
      discovery_panel == "Lung antiviral IFN module" ~ 1,
      discovery_panel == "Tonsil complement/coagulation module" ~ 2,
      TRUE ~ 99
    ),
    gene_label = gene,
    candidate_class = case_when(
      manuscript_module == "Conserved antiviral/IFN core" ~ "Conserved antiviral/IFN core",
      manuscript_module == "Mainly in vivo immune/progression module" ~ "In vivo IFN/progression support",
      manuscript_module == "In vivo complement/coagulation disease module" ~ "Complement/coagulation disease module",
      TRUE ~ manuscript_module
    )
  ) %>%
  arrange(panel_order, desc(priority_score), desc(abs_kME)) %>%
  select(
    discovery_panel, tissue, module_label, module_color, module_numeric,
    gene, candidate_class, category, manuscript_module, evidence_strength,
    kME, abs_kME, infected_module_cor, infected_module_p,
    dpi_module_cor, dpi_module_p, priority_score, n_significant_up,
    n_huvec_significant_up, n_in_vivo_significant_up,
    conserved_call, wgcna_support_class, discovery_interpretation
  )

if (nrow(focused) == 0) {
  stop("No focused candidates were selected. Check module numbers/colors and 03c thresholds.")
}

module_summary <- focused %>%
  distinct(
    discovery_panel, tissue, module_label, module_color, module_numeric,
    infected_module_cor, infected_module_p, dpi_module_cor, dpi_module_p,
    discovery_interpretation
  ) %>%
  mutate(
    n_focused_genes = map_int(discovery_panel, ~sum(focused$discovery_panel == .x)),
    manuscript_sentence = case_when(
      tissue == "Lung" ~ paste0(
        "The lung tan module (ME12) was strongly associated with infection (r = ",
        round(infected_module_cor, 2), ") and DPI (r = ", round(dpi_module_cor, 2),
        "), and contained high-membership antiviral/IFN candidates."
      ),
      tissue == "Tonsil" ~ paste0(
        "The tonsil blue module (ME2) was strongly associated with infection (r = ",
        round(infected_module_cor, 2), ") and DPI (r = ", round(dpi_module_cor, 2),
        "), and contained high-membership complement/coagulation candidates."
      ),
      TRUE ~ discovery_interpretation
    )
  )

readr::write_csv(
  focused,
  file.path(table_dir, "focused_wgcna_candidate_modules.csv")
)

wb <- createWorkbook()
addWorksheet(wb, "focused_candidates")
writeData(wb, "focused_candidates", focused)
freezePane(wb, "focused_candidates", firstRow = TRUE)
setColWidths(wb, "focused_candidates", cols = 1:ncol(focused), widths = "auto")

addWorksheet(wb, "module_summary")
writeData(wb, "module_summary", module_summary)
freezePane(wb, "module_summary", firstRow = TRUE)
setColWidths(wb, "module_summary", cols = 1:ncol(module_summary), widths = "auto")

saveWorkbook(
  wb,
  file.path(table_dir, "focused_wgcna_candidate_modules.xlsx"),
  overwrite = TRUE
)

###############################################################################
# Publication-style focused panel
###############################################################################

panel_palette <- c(
  "Conserved antiviral/IFN core" = "#2C7FB8",
  "In vivo IFN/progression support" = "#7B3294",
  "Complement/coagulation disease module" = "#D95F02"
)

make_module_plot <- function(df, panel_label) {
  plot_df <- df %>%
    arrange(abs_kME) %>%
    mutate(gene_label = factor(gene, levels = unique(gene)))  # fix: gene_label is dropped by the select() above; build it from gene

  module_r <- plot_df %>%
    summarise(
      infected_r = dplyr::first(infected_module_cor),
      dpi_r = dplyr::first(dpi_module_cor)
    )

  ggplot(plot_df, aes(x = abs_kME, y = gene_label)) +
    geom_segment(
      aes(x = 0.78, xend = abs_kME, yend = gene_label),
      linewidth = 0.45,
      color = "grey72"
    ) +
    geom_point(
      aes(fill = candidate_class, size = priority_score),
      shape = 21,
      color = "grey15",
      linewidth = 0.35,
      alpha = 0.96
    ) +
    geom_vline(
      xintercept = 0.8,
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey35"
    ) +
    annotate(
      "label",
      x = 0.807,
      y = Inf,
      hjust = 0,
      vjust = 1.2,
      label = paste0(
        "module-infection r = ", round(module_r$infected_r, 2),
        "\nmodule-DPI r = ", round(module_r$dpi_r, 2)
      ),
      size = 3.2,
      label.size = 0.25,
      color = "grey15",
      fill = "white"
    ) +
    scale_x_continuous(
      limits = c(0.78, 1.005),
      breaks = c(0.80, 0.85, 0.90, 0.95, 1.00),
      labels = number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    scale_fill_manual(values = panel_palette, drop = FALSE) +
    scale_size_continuous(range = c(3.2, 7.2), guide = "none") +
    labs(
      title = panel_label,
      x = "|kME| module membership",
      y = NULL,
      fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15, color = "grey10"),
      axis.text.y = element_text(face = "plain", color = "grey18", size = 10.5),
      axis.text.x = element_text(color = "grey20", size = 10.5),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      plot.margin = margin(16, 18, 10, 10)
    )
}

p_lung <- make_module_plot(
  focused %>% filter(discovery_panel == "Lung antiviral IFN module"),
  "A"
)

p_tonsil <- make_module_plot(
  focused %>% filter(discovery_panel == "Tonsil complement/coagulation module"),
  "B"
)

combined_panel <- p_lung / p_tonsil +
  plot_layout(guides = "collect", heights = c(1.35, 0.9)) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 4)
  ) &
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

png(
  file.path(figure_dir, "focused_wgcna_candidate_modules_panel.png"),
  width = 3600, height = 3200, res = 300
)
print(combined_panel)
dev.off()

###############################################################################
# Compact heatmap-style summary
###############################################################################

heatmap_df <- focused %>%
  mutate(
    label = paste(tissue, gene, sep = ": "),
    module_infection_r = infected_module_cor,
    module_DPI_r = dpi_module_cor
  ) %>%
  select(label, abs_kME, module_infection_r, module_DPI_r) %>%
  pivot_longer(
    cols = c(abs_kME, module_infection_r, module_DPI_r),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("abs_kME", "module_infection_r", "module_DPI_r"),
      labels = c("|kME|", "Module-infection r", "Module-DPI r")
    ),
    label = factor(label, levels = rev(unique(focused %>% mutate(label = paste(tissue, gene, sep = ": ")) %>% pull(label))))
  )

png(
  file.path(figure_dir, "focused_wgcna_candidate_modules_heatmap.png"),
  width = 2600, height = 2200, res = 300
)
ggplot(heatmap_df, aes(x = metric, y = label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.8, color = "grey10") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Value"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", color = "grey15"),
    axis.text.y = element_text(face = "plain", color = "grey20", size = 10),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(8, 12, 8, 8)
  )
dev.off()

message("Focused WGCNA discovery figure complete.")
message("Table: ", file.path(table_dir, "focused_wgcna_candidate_modules.csv"))
message("Panel: ", file.path(figure_dir, "focused_wgcna_candidate_modules_panel.png"))
message("Heatmap: ", file.path(figure_dir, "focused_wgcna_candidate_modules_heatmap.png"))
