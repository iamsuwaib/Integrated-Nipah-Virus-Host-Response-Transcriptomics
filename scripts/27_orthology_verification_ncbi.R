###############################################################################
# Orthology verification, second curated resource: NCBI Gene orthologs (RefSeq)
#
# Script 26 queried Ensembl Compara. The GSE310471 count matrix was annotated
# with RefSeq gene symbols (Chlorocebus sabaeus), and Ensembl's vervet
# annotation is less complete than RefSeq's, so some candidate genes had no
# Ensembl orthologue entry. This script therefore repeats the check against
# NCBI Gene's curated ortholog assignments, which are built on the same RefSeq
# annotation as the dataset, and combines both resources.
#
# Resources (downloaded once, cached):
#   gene_orthologs.gz                       NCBI curated ortholog pairs
#   Homo_sapiens.gene_info.gz               human gene symbols / IDs
# and, through NCBI E-utilities (small queries, no bulk download), the symbols of
# the AGM ortholog genes (taxid 60711) and whether an AGM gene with the
# identical symbol exists in RefSeq.
#
# For each candidate gene (same set as script 26) the script reports
#   - the human NCBI GeneID,
#   - the AGM ortholog GeneID(s) and symbol(s) in NCBI,
#   - the ortholog multiplicity in both directions (one-to-one when exactly one
#     AGM gene is assigned to the human gene AND exactly one human gene to that
#     AGM gene),
#   - whether the AGM ortholog symbol equals the human symbol used for matching,
#   - whether an AGM gene with the identical symbol exists in RefSeq at all,
# and combines this with the Ensembl result from script 26.
#
# Outputs (advanced_analyses/tables):
#   orthology_verification_ncbi.csv
#   orthology_verification_combined.csv
#   orthology_verification_combined_summary.csv
###############################################################################

options(stringsAsFactors = FALSE, timeout = 3600)
for (pkg in c("tidyverse", "httr", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
suppressPackageStartupMessages({ library(tidyverse); library(httr); library(jsonlite) })

project_dir <- "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
table_dir <- file.path(project_dir, "advanced_analyses", "tables")
cache_dir <- file.path(project_dir, "advanced_analyses", "external_resources")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

urls <- c(
  gene_orthologs = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
  human_info = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
)
files <- file.path(cache_dir, basename(urls)); names(files) <- names(urls)
for (n in names(urls)) {
  if (!file.exists(files[[n]])) {
    message("Downloading ", urls[[n]])
    download.file(urls[[n]], files[[n]], mode = "wb", quiet = FALSE)
  }
}

HUMAN <- 9606L; AGM <- 60711L

hinfo <- readr::read_tsv(files[["human_info"]], comment = "", show_col_types = FALSE,
                         col_select = c(GeneID, Symbol, Synonyms, type_of_gene)) %>%
  filter(type_of_gene == "protein-coding" | TRUE)

orth <- readr::read_tsv(files[["gene_orthologs"]], comment = "", show_col_types = FALSE,
                        col_names = c("tax_id", "GeneID", "Relationship", "Other_tax_id", "Other_GeneID"),
                        skip = 1,
                        col_types = cols(.default = col_integer(), Relationship = col_character())) %>%
  filter((tax_id == HUMAN & Other_tax_id == AGM))
message("Human-AGM ortholog pairs in NCBI: ", nrow(orth))

n_agm_per_human <- orth %>% count(GeneID, name = "n_agm_for_human")
n_human_per_agm <- orth %>% count(Other_GeneID, name = "n_human_for_agm")
orth <- orth %>% left_join(n_agm_per_human, by = "GeneID") %>% left_join(n_human_per_agm, by = "Other_GeneID")

# candidate set (identical to script 26)
short1 <- readr::read_csv(file.path(table_dir, "candidate_biomarker_shortlist.csv"), show_col_types = FALSE)
sec31  <- readr::read_csv(file.path(table_dir, "secretome_surfaceome_candidate_table_refined.csv"), show_col_types = FALSE)
wgcna  <- readr::read_csv(file.path(table_dir, "focused_wgcna_candidate_modules.csv"), show_col_types = FALSE)
cand <- tibble(gene = sort(unique(c(short1$gene, sec31$gene, wgcna$gene))))

hsym <- hinfo %>% filter(Symbol %in% cand$gene) %>% select(human_geneid = GeneID, gene = Symbol)

# --- NCBI E-utilities helpers (<= 3 requests/second without an API key) ---
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
eget <- function(endpoint, query, tries = 4) {
  for (k in seq_len(tries)) {
    Sys.sleep(0.4)
    r <- tryCatch(GET(paste0(eutils, endpoint), query = c(query, retmode = "json"), timeout(60)), error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200)
      return(fromJSON(content(r, as = "text", encoding = "UTF-8"), simplifyVector = FALSE))
    Sys.sleep(1.5 * k)
  }
  NULL
}
# AGM gene with the identical official symbol present in RefSeq?
agm_same_symbol <- map_dfr(cand$gene, function(g) {
  j <- eget("esearch.fcgi", list(db = "gene", term = paste0(g, "[Gene Name] AND 60711[Taxonomy ID]")))
  ids <- if (is.null(j)) character() else unlist(j$esearchresult$idlist)
  tibble(gene = g, agm_geneid_same_symbol = if (length(ids)) as.integer(ids[1]) else NA_integer_)
})
# symbols of AGM ortholog genes
agm_ids <- unique(orth$Other_GeneID[orth$GeneID %in% hsym$human_geneid])
ainfo <- map_dfr(split(agm_ids, ceiling(seq_along(agm_ids) / 100)), function(ch) {
  j <- eget("esummary.fcgi", list(db = "gene", id = paste(ch, collapse = ",")))
  if (is.null(j)) return(tibble(GeneID = integer(), Symbol = character()))
  res <- j$result; res$uids <- NULL
  tibble(GeneID = as.integer(names(res)), Symbol = map_chr(res, ~ .x$name %||% NA_character_))
})

ncbi <- cand %>%
  left_join(hsym, by = "gene") %>%
  left_join(orth %>% select(human_geneid = GeneID, agm_geneid = Other_GeneID, n_agm_for_human, n_human_for_agm),
            by = "human_geneid", relationship = "many-to-many") %>%
  left_join(ainfo %>% select(agm_geneid = GeneID, agm_symbol = Symbol), by = "agm_geneid") %>%
  left_join(agm_same_symbol, by = "gene") %>%
  mutate(
    ncbi_status = case_when(
      is.na(human_geneid) ~ "human symbol not found in NCBI",
      is.na(agm_geneid) & !is.na(agm_geneid_same_symbol) ~ "no NCBI ortholog pair, but AGM gene with identical symbol exists",
      is.na(agm_geneid) ~ "no NCBI ortholog found",
      n_agm_for_human == 1 & n_human_for_agm == 1 & toupper(agm_symbol) == toupper(gene) ~ "one-to-one, symbol matches",
      n_agm_for_human == 1 & n_human_for_agm == 1 ~ "one-to-one, AGM symbol differs",
      TRUE ~ "one-to-many or many-to-many"
    )
  ) %>%
  select(gene, human_geneid, agm_geneid, agm_symbol, n_agm_for_human, n_human_for_agm,
         agm_geneid_same_symbol, ncbi_status)
readr::write_csv(ncbi, file.path(table_dir, "orthology_verification_ncbi.csv"))

# combine with Ensembl (script 26)
ens <- readr::read_csv(file.path(table_dir, "orthology_verification_candidates.csv"), show_col_types = FALSE) %>%
  group_by(gene) %>%
  summarise(ensembl_status = paste(unique(status), collapse = "; "),
            ensembl_agm_symbol = paste(unique(na.omit(agm_symbol)), collapse = ","), .groups = "drop")

combined <- ncbi %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.x)), collapse = ",")), .groups = "drop") %>%
  left_join(ens, by = "gene") %>%
  mutate(
    ncbi_ok = ncbi_status == "one-to-one, symbol matches",
    ensembl_ok = ensembl_status == "one-to-one, symbol matches",
    overall = case_when(
      ncbi_ok & ensembl_ok ~ "confirmed by both resources",
      ncbi_ok ~ "confirmed by NCBI (RefSeq) only",
      ensembl_ok ~ "confirmed by Ensembl only",
      TRUE ~ "not confirmed as one-to-one by either resource"
    )
  )
readr::write_csv(combined, file.path(table_dir, "orthology_verification_combined.csv"))

summ <- combined %>% count(overall, name = "n_genes")
readr::write_csv(summ, file.path(table_dir, "orthology_verification_combined_summary.csv"))
print(summ)
message("Genes not confirmed by both resources:")
print(combined %>% filter(overall != "confirmed by both resources") %>%
        select(gene, ncbi_status, agm_symbol, ensembl_status, overall), n = Inf, width = Inf)
message("Done")
