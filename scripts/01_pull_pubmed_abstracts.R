#!/usr/bin/env Rscript

# 01_pull_pubmed_abstracts.R
#
# Retrieves PubMed records for searches combining archaic introgression terms
# and neurobehavioural trait terms (2010 onwards), then downloads XML records
# and extracts metadata + abstracts.
#
# Outputs (written under ./data/):
#   data/logs/pubmed_query_log.csv
#   data/logs/pubmed_pmids_all.csv
#   data/raw/pubmed_abstracts_raw.csv
#   data/processed/pubmed_with_abstracts.csv
#   data/processed/pubmed_no_abstracts.csv
#
# Author: David Pulford
# Year: 2026

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(rentrez)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(xml2)
  library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

safe_text <- function(node) {
  if (inherits(node, "xml_missing") || length(node) == 0) return(NA_character_)
  txt <- xml_text(node)
  if (is.null(txt) || txt == "") return(NA_character_)
  txt
}

# ---- Paths ----
BASE_DIR <- file.path(getwd(), "data")

dir.create(file.path(BASE_DIR, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_DIR, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_DIR, "processed"), recursive = TRUE, showWarnings = FALSE)

OUT_RAW   <- file.path(BASE_DIR, "raw", "pubmed_abstracts_raw.csv")
OUT_LOG   <- file.path(BASE_DIR, "logs", "pubmed_query_log.csv")
OUT_PMIDS <- file.path(BASE_DIR, "logs", "pubmed_pmids_all.csv")
OUT_WITH  <- file.path(BASE_DIR, "processed", "pubmed_with_abstracts.csv")
OUT_NO    <- file.path(BASE_DIR, "processed", "pubmed_no_abstracts.csv")

# ---- Query terms ----
TRAIT_TERMS <- c(
  "cognition", "cognitive", "memory", "executive function", "language",
  "psychiatric", "schizophrenia", "bipolar", "depression", "anxiety",
  "autism", "ADHD", "behaviour", "behavior", "risk taking", "impulsivity",
  "addiction", "smoking", "nicotine", "circadian", "chronotype", "sleep",
  "pain", "sociability", "social behaviour", "social behavior"
)

ARCHAIC_TERMS <- c(
  "\"Neanderthal introgression\"",
  "\"archaic introgression\"",
  "introgressed",
  "archaic haplotype",
  "Neanderthal-derived",
  "Neandertal",
  "Denisovan"
)

DATE_FILTER <- "\"2010/01/01\"[Date - Publication] : \"3000\"[Date - Publication]"
archaic_block <- paste0("(", paste(ARCHAIC_TERMS, collapse = " OR "), ")")

build_query <- function(trait) {
  trait_q <- if (str_detect(trait, "\\s")) paste0("\"", trait, "\"") else trait
  paste0(archaic_block, " AND ", trait_q, " AND ", DATE_FILTER)
}

queries <- tibble(
  trait = TRAIT_TERMS,
  query = vapply(TRAIT_TERMS, build_query, FUN.VALUE = character(1))
)

# ---- PubMed search ----
search_pubmed <- function(q, retmax = 5000) {
  res <- tryCatch(
    entrez_search(db = "pubmed", term = q, retmax = retmax),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(count = NA_integer_, ids = character(0)))
  list(count = res$count, ids = res$ids)
}

# ---- XML fetch + parse ----
fetch_xml_records <- function(pmids) {
  if (length(pmids) == 0) return(tibble())

  xml_txt <- tryCatch(
    entrez_fetch(db = "pubmed", id = pmids, rettype = "xml", retmode = "text"),
    error = function(e) NULL
  )
  if (is.null(xml_txt)) return(tibble())

  doc <- read_xml(xml_txt)
  articles <- xml_find_all(doc, ".//PubmedArticle")

  map_dfr(articles, function(article) {
    pmid <- safe_text(xml_find_first(article, ".//MedlineCitation/PMID"))

    title <- safe_text(xml_find_first(article, ".//ArticleTitle"))
    journal <- safe_text(xml_find_first(article, ".//Journal/Title"))

    year <- safe_text(xml_find_first(article, ".//JournalIssue/PubDate/Year"))
    if (is.na(year)) {
      medline <- safe_text(xml_find_first(article, ".//JournalIssue/PubDate/MedlineDate"))
      year <- str_extract(medline %||% "", "\\d{4}")
      if (identical(year, character(0)) || year == "") year <- NA_character_
    }

    auth_nodes <- xml_find_all(article, ".//AuthorList/Author")
    authors <- if (length(auth_nodes) == 0) {
      NA_character_
    } else {
      auth_names <- map_chr(auth_nodes, function(a) {
        ln <- safe_text(xml_find_first(a, ".//LastName"))
        inits <- safe_text(xml_find_first(a, ".//Initials"))
        nm <- str_squish(paste(ln, inits))
        if (nm == "") NA_character_ else nm
      })
      out <- paste(na.omit(auth_names), collapse = "; ")
      if (out == "") NA_character_ else out
    }

    abs_nodes <- xml_find_all(article, ".//Abstract/AbstractText")
    abstract <- if (length(abs_nodes) == 0) NA_character_ else paste(xml_text(abs_nodes), collapse = " ")

    tibble(
      pmid = pmid,
      title = title,
      journal = journal,
      year = year,
      authors = authors,
      abstract = abstract
    )
  })
}

# ---- Run searches ----
message("PubMed retrieval started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

log_rows <- list()
all_ids <- character(0)

for (i in seq_len(nrow(queries))) {
  tr <- queries$trait[i]
  q <- queries$query[i]

  message(sprintf("[%d/%d] %s", i, nrow(queries), tr))
  out <- search_pubmed(q)

  log_rows[[i]] <- tibble(
    trait = tr,
    query = q,
    count = out$count
  )

  if (length(out$ids) > 0) all_ids <- c(all_ids, out$ids)
  Sys.sleep(0.34)
}

log_df <- bind_rows(log_rows)
write_csv(log_df, OUT_LOG)

all_ids <- unique(all_ids)
write_csv(tibble(pmid = all_ids), OUT_PMIDS)

message("Unique PMIDs: ", length(all_ids))

# ---- Fetch records ----
BATCH_SIZE <- 100
rec_list <- list()

for (k in seq(1, length(all_ids), by = BATCH_SIZE)) {
  batch <- all_ids[k:min(k + BATCH_SIZE - 1, length(all_ids))]
  rec_list[[length(rec_list) + 1]] <- fetch_xml_records(batch)
  Sys.sleep(0.34)
}

final_df <- bind_rows(rec_list) %>%
  distinct(pmid, .keep_all = TRUE)

write_csv(final_df, OUT_RAW)

with_abs <- final_df %>% filter(!is.na(abstract) & str_squish(abstract) != "")
no_abs <- final_df %>% filter(is.na(abstract) | str_squish(abstract) == "")

write_csv(with_abs, OUT_WITH)
write_csv(no_abs, OUT_NO)

message("Saved: ", OUT_RAW)
message("With abstracts: ", nrow(with_abs))
message("No abstracts: ", nrow(no_abs))
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

