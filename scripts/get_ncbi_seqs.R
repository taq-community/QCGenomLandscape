#!/usr/bin/env Rscript
# Script de test rapide pour vérifier les APIs
# Teste avec quelques espèces représentatives
library(tidyverse)
library(httr)
library(jsonlite)
library(logger)
library(rentrez)

# Configure logger to write only to file (no console output)
log_dir <- "logs"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

log_file <- file.path(log_dir, sprintf("ncbi_queries_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_appender(appender_file(log_file))
log_threshold(INFO)

log_info("Logger initialized. Writing to: {log_file}")

query_primers <- read.csv2("data/primers_map_group_bdqc_list_01122025.csv")
  
qc_species <- read.csv("data/bdqc_list_01122025.csv") |>
  filter(rank == "species") |>
  mutate(group = tolower(group_en)) |>
  left_join(query_primers |> mutate(group = tolower(group)), by = c("group")) |>
  mutate(query = ifelse(!is.na(query_marker), glue::glue("{species}[Organism] AND {query_marker} AND voucher[Title]"), NA))

queries <- qc_species |>
  pull(query) |>
  na.omit() |>
  unique()

# Initialize logger
log_info("Starting NCBI queries for {length(queries)} species")
set_entrez_key(Sys.getenv("NCBI_API_KEY"))

deficient_queries <- list()
high_id_queries <- list()

ncbi_results <- map_df(seq_along(queries), \(i) {
  q <- queries[14010]
  log_info("Query {i}/{length(queries)}: {q}")

  # Wrap entire query process in tryCatch
  tryCatch(
    {
      # First search to get total count - wrapped in tryCatch
      results <- tryCatch(
        rentrez::entrez_search(db = "nucleotide", term = q, retmax = 0),
        error = function(e) {
          log_error("Error in initial search for query {i}: {e$message}")
          deficient_queries[[length(deficient_queries) + 1]] <<- list(
            query_index = i,
            query = q,
            error_type = "entrez_search_count",
            error_message = e$message,
            timestamp = Sys.time()
          )
          return(NULL)
        }
      )

      if (is.null(results)) {
        return(tibble())
      }

      total_count <- results$count

      if (total_count == 0) {
        log_warn("No results found for query")
        return(tibble())
      }

      log_info("Found {total_count} total sequences, retrieving all IDs...")

      # Get all IDs using retmax - wrapped in tryCatch
      results <- tryCatch(
        rentrez::entrez_search(db = "nucleotide", term = q, retmax = 5000),
        error = function(e) {
          log_error("Error retrieving IDs for query {i}: {e$message}")
          deficient_queries[[length(deficient_queries) + 1]] <<- list(
            query_index = i,
            query = q,
            error_type = "entrez_search_ids",
            error_message = e$message,
            timestamp = Sys.time()
          )
          return(NULL)
        }
      )

      retrieved_ids <- length(results$id)

      log_success("Retrieved {retrieved_ids} IDs")

      # Store high ID queries (> 500)
      if (retrieved_ids > 500) {
        high_id_queries[[length(high_id_queries) + 1]] <<- list(
          query_index = i,
          query = q,
          id_count = retrieved_ids,
          timestamp = Sys.time()
        )
        log_warn("High ID count ({retrieved_ids} > 500) - query stored")
      }

      # Batch summaries to avoid HTTP 414 error (max ~200 IDs per request)
      batch_size <- 200
      all_summaries <- list()

      for (batch_start in seq(1, length(results$id), by = batch_size)) {
        batch_end <- min(batch_start + batch_size - 1, length(results$id))
        batch_ids <- results$id[batch_start:batch_end]

        log_info("Fetching summaries {batch_start}-{batch_end} of {length(results$id)}")

        # Wrap summary fetch in tryCatch
        batch_summary <- tryCatch(
          {
            rentrez::entrez_summary(db = "nucleotide", id = batch_ids)
          },
          error = function(e) {
            log_error("Error fetching summaries for batch {batch_start}-{batch_end}: {e$message}")
            deficient_queries[[length(deficient_queries) + 1]] <<- list(
              query_index = i,
              query = q,
              error_type = "entrez_summary",
              error_message = e$message,
              batch_range = paste0(batch_start, "-", batch_end),
              batch_ids = batch_ids,
              timestamp = Sys.time()
            )
            NULL
          }
        )

        # Store batch results if successful
        if (!is.null(batch_summary)) {
          if (length(batch_ids) == 1) {
            all_summaries[[length(all_summaries) + 1]] <- batch_summary
          } else {
            all_summaries <- c(all_summaries, batch_summary)
          }
        }
      }

      results_df <- map_df(all_summaries, function(x) {
        if (is.list(x)) {
          # Split subtype and subname - wrap in tryCatch for atomic vectors
          subtypes <- tryCatch(
            {
              strsplit(x$subtype, "\\|")[[1]]
            },
            error = function(e) character(0)
          )

          subnames <- tryCatch(
            {
              strsplit(x$subname, "\\|")[[1]]
            },
            error = function(e) character(0)
          )

          # Create named list for subtypes
          subtype_values <- if (length(subtypes) > 0 && length(subnames) > 0) {
            setNames(as.list(subnames), subtypes)
          } else {
            list()
          }

          tibble(
            uid = x$uid %||% NA,
            accession = x$accessionversion %||% NA,
            title = x$title %||% NA,
            taxid = x$taxid %||% NA,
            organism = x$organism %||% NA,
            moltype = x$moltype %||% NA,
            topology = x$topology %||% NA,
            genome = x$genome %||% NA,
            slen = x$slen %||% NA,
            createdate = x$createdate %||% NA,
            updatedate = x$updatedate %||% NA,
            specimen_voucher = subtype_values$specimen_voucher %||% NA,
            country = subtype_values$country %||% NA,
            lat_lon = subtype_values$lat_lon %||% NA,
            collection_date = subtype_values$collection_date %||% NA
          )
        }
      }) |> mutate(query = q)
    },
    error = function(e) {
      log_error("Error processing query {i}: {e$message}")
      deficient_queries[[length(deficient_queries) + 1]] <<- list(
        query_index = i,
        query = q,
        error_type = "entrez_search",
        error_message = e$message,
        timestamp = Sys.time()
      )
      return(tibble())
    }
  )
}, .progress = TRUE) |> filter(!if_all(everything(), is.na)) 

saveRDS(ncbi_results, "results/ncbi_results.rds")
saveRDS(deficient_queries, "results/deficient_queries.rds")
saveRDS(high_id_queries, "results/high_id_queries.rds")

log_success("Completed! Retrieved {nrow(ncbi_results)} total sequences")
