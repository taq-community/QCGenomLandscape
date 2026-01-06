#!/usr/bin/env Rscript
# Script to retrieve sequences from BOLD Systems for QC species
library(tidyverse)
library(bold)
library(logger)

# Configure logger to write only to file (no console output)
log_dir <- "logs"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

log_file <- file.path(log_dir, sprintf("bold_queries_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_appender(appender_file(log_file))
log_threshold(INFO)

log_info("Logger initialized. Writing to: {log_file}")

# Load species list
qc_species <- read.csv("data/bdqc_list_01122025.csv") |>
  filter(rank == "species")

# Initialize logger
log_info("Starting BOLD queries for {nrow(qc_species)} species")

deficient_queries <- list()

bold_results <- map_df(seq_len(nrow(qc_species))[100:120], \(i) {
    i = 102
    species_name <- qc_species$species[i]
    log_info("Query {i}/{nrow(qc_species)}: {species_name}")

    # Wrap entire query process in tryCatch
    tryCatch(
        {
            # Query BOLD for sequences - wrapped in tryCatch
            seq_data <- tryCatch(
                bold::bold_specimens(taxon = species_name, format = "tsv", cleanData = TRUE),
                error = function(e) {
                    log_error("Error in BOLD query for {species_name}: {e$message}")
                    deficient_queries[[length(deficient_queries) + 1]] <<- list(
                        query_index = i,
                        species = species_name,
                        error_type = "bold_seqspec",
                        error_message = e$message,
                        timestamp = Sys.time()
                    )
                    return(NULL)
                }
            )

            if (is.null(seq_data) || nrow(seq_data) == 0) {
                log_warn("No results found for {species_name}")
                return(tibble())
            }

            # Add species name to results and select relevant columns
            seq_data <- seq_data |> mutate(query_species = species_name)
            return(seq_data)
        },
        error = function(e) {
            log_error("Error processing query {i} ({species_name}): {e$message}")
            deficient_queries[[length(deficient_queries) + 1]] <<- list(
                query_index = i,
                species = species_name,
                error_type = "general",
                error_message = e$message,
                timestamp = Sys.time()
            )
            return(tibble())
        }
    )

    # Add small delay to avoid overwhelming BOLD API
    Sys.sleep(0.5)
}, .progress = TRUE)

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

saveRDS(bold_results, "results/bold_results.rds")
saveRDS(deficient_queries, "results/deficient_queries.rds")

log_success("Completed! Retrieved {nrow(bold_results)} total sequences from {length(unique(bold_results$query_species))} species")
