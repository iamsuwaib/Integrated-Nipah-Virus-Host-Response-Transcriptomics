###############################################################################
# One-to-one orthology verification of the final candidate gene set
# (human <-> African green monkey, Chlorocebus sabaeus)
#
# The cross-dataset integration matches human (HUVEC) and African green monkey
# (GSE310471) results by gene symbol. To verify that this symbol matching is
# valid for the genes the manuscript actually relies on, every gene in the
# final candidate set is checked against a curated orthology resource, Ensembl
# Compara (REST API, rest.ensembl.org), which classifies each human gene's
# Chlorocebus sabaeus orthologue(s) as one-to-one, one-to-many or many-to-many
# and gives sequence identity and a gene-order/whole-genome-alignment based
# high-confidence flag where available.
#
# Final candidate set = union of
#   (i)   the Table 1 candidate shortlist           (candidate_biomarker_shortlist.csv)
#   (ii)  the 31-member secretome/surfaceome catalogue
#                                                   (secretome_surfaceome_candidate_table_refined.csv)
#   (iii) the genes displayed in the focused WGCNA modules (Figure 4A)
#                                                   (focused_wgcna_candidate_modules.csv)
#
# For each gene the output reports: Ensembl human gene ID, AGM orthologue
# Ensembl ID and symbol, orthology type, percent identity (both directions),
# the high-confidence flag if provided, whether the AGM symbol equals the human
# symbol used for matching, and an overall status. Ambiguous mappings
# (one-to-many, many-to-many, no orthologue, or symbol mismatch) are listed
# explicitly.
#
# Requires internet access (Ensembl REST allows ~15 requests/second; this
# script sleeps between calls). Runtime: about 2-4 minutes.
#
# Outputs (advanced_analyses/tables):
#   orthology_verification_candidates.csv
#   orthology_verification_summary.csv
###############################################################################

options(stringsAsFactors = FALSE)
for (pkg in c("tidyverse", "httr", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
})

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
table_dir <- file.path(project_dir, "advanced_analyses", "tables")

short1 <- readr::read_csv(file.path(table_dir, "candidate_biomarker_shortlist.csv"), show_col_types = FALSE)
sec31  <- readr::read_csv(file.path(table_dir, "secretome_surfaceome_candidate_table_refined.csv"), show_col_types = FALSE)
wgcna  <- readr::read_csv(file.path(table_dir, "focused_wgcna_candidate_modules.csv"), show_col_types = FALSE)

cand <- tibble(gene = sort(unique(c(short1$gene, sec31$gene, wgcna$gene)))) %>%
  mutate(
    in_table1_shortlist     = gene %in% short1$gene,
    in_secretome_catalogue  = gene %in% sec31$gene,
    in_wgcna_display        = gene %in% wgcna$gene
  )
message("Final candidate set: ", nrow(cand), " genes")

SERVER <- "https://rest.ensembl.org"

get_json <- function(path, tries = 5) {
  for (k in seq_len(tries)) {
    resp <- tryCatch(
      GET(paste0(SERVER, path), accept_json(), timeout(60)),
      error = function(e) NULL
    )
    if (!is.null(resp)) {
      if (status_code(resp) == 200) return(fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE))
      if (status_code(resp) == 429) {
        wait <- suppressWarnings(as.numeric(headers(resp)[["retry-after"]]))
        Sys.sleep(ifelse(is.na(wait), 2, wait + 0.5)); next
      }
      if (status_code(resp) %in% c(400, 404)) return(NULL)  # gene not found / no homology
    }
    Sys.sleep(1.5 * k)
  }
  NULL
}

# AGM gene symbol lookup (batch via POST)
lookup_symbols <- function(ids) {
  ids <- unique(ids[!is.na(ids)])
  if (!length(ids)) return(tibble(agm_id = character(), agm_symbol = character()))
  resp <- POST(paste0(SERVER, "/lookup/id"), accept_json(), content_type_json(),
               body = toJSON(list(ids = as.list(ids)), auto_unbox = TRUE), timeout(120))
  if (status_code(resp) != 200) return(tibble(agm_id = ids, agm_symbol = NA_character_))
  res <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  tibble(agm_id = names(res),
         agm_symbol = map_chr(res, ~ if (is.null(.x$display_name)) NA_character_ else .x$display_name))
}

query_gene <- function(sym) {
  Sys.sleep(0.15)
  j <- get_json(paste0("/homology/symbol/homo_sapiens/", sym,
                       "?target_species=chlorocebus_sabaeus;type=orthologues;format=full"))
  if (is.null(j) || length(j$data) == 0 || length(j$data[[1]]$homologies) == 0) {
    return(tibble(gene = sym, human_id = NA_character_, agm_id = NA_character_,
                  orthology_type = "no_orthologue_found", perc_id_agm_to_human = NA_real_,
                  perc_id_human_to_agm = NA_real_, high_confidence = NA, n_agm_orthologues = 0L))
  }
  h <- j$data[[1]]$homologies
  n <- length(h)
  map_dfr(h, function(x) {
    tibble(gene = sym,
           human_id = x$source$id %||% j$data[[1]]$id,
           agm_id = x$target$id,
           orthology_type = x$type,
           perc_id_agm_to_human = as.numeric(x$target$perc_id %||% NA),
           perc_id_human_to_agm = as.numeric(x$source$perc_id %||% NA),
           high_confidence = if (is.null(x$is_high_confidence)) NA else as.logical(as.integer(x$is_high_confidence)),
           n_agm_orthologues = n)
  })
}

raw <- map_dfr(cand$gene, function(g) { message(g); query_gene(g) })
syms <- lookup_symbols(raw$agm_id)

res <- raw %>%
  left_join(syms, by = "agm_id") %>%
  left_join(cand, by = "gene") %>%
  mutate(
    symbol_match = !is.na(agm_symbol) & toupper(agm_symbol) == toupper(gene),
    status = case_when(
      orthology_type == "no_orthologue_found" ~ "no orthologue found",
      orthology_type == "ortholog_one2one" & symbol_match ~ "one-to-one, symbol matches",
      orthology_type == "ortholog_one2one" & !symbol_match ~ "one-to-one, AGM symbol differs or is unnamed",
      orthology_type == "ortholog_one2many" ~ "one-to-many (ambiguous)",
      orthology_type == "ortholog_many2many" ~ "many-to-many (ambiguous)",
      TRUE ~ paste("other:", orthology_type)
    )
  ) %>%
  select(gene, in_table1_shortlist, in_secretome_catalogue, in_wgcna_display,
         human_id, agm_id, agm_symbol, orthology_type, n_agm_orthologues,
         perc_id_agm_to_human, perc_id_human_to_agm, high_confidence, symbol_match, status)

readr::write_csv(res, file.path(table_dir, "orthology_verification_candidates.csv"))

gene_level <- res %>%
  group_by(gene) %>%
  summarise(clean = any(status == "one-to-one, symbol matches") && n() == 1,
            status = paste(unique(status), collapse = "; "), .groups = "drop")

summary_tbl <- tibble(
  metric = c("candidate genes checked",
             "one-to-one orthologue with matching symbol",
             "one-to-one orthologue, symbol differs or unnamed",
             "one-to-many or many-to-many",
             "no orthologue found"),
  n = c(nrow(gene_level),
        sum(gene_level$status == "one-to-one, symbol matches"),
        sum(grepl("symbol differs", gene_level$status)),
        sum(grepl("many", gene_level$status)),
        sum(grepl("no orthologue", gene_level$status)))
)
readr::write_csv(summary_tbl, file.path(table_dir, "orthology_verification_summary.csv"))
print(summary_tbl)
amb <- res %>% filter(status != "one-to-one, symbol matches")
if (nrow(amb)) { message("Ambiguous or unresolved mappings:"); print(amb %>% select(gene, agm_id, agm_symbol, orthology_type, status), n = Inf) }
message("Done")
