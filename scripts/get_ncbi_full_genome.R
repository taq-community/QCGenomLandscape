#!/usr/bin/env Rscript
# Query NCBI for full genomes (nuclear or mitochondrial) for each species
library(tidyverse)
library(rentrez)
library(logger)

# Configure logger
log_dir <- "logs"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

log_file <- file.path(log_dir, sprintf("ncbi_full_genome_queries_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_appender(appender_file(log_file))
log_threshold(INFO)

log_info("Logger initialized. Writing to: {log_file}")

entrez_dbs()

# Set NCBI API key
set_entrez_key(Sys.getenv("NCBI_API_KEY"))

# Load species list
qc_species <- read.csv("data/bdqc_list_01122025.csv") |>
  filter(rank == "species") |>
  distinct(species) |>
  pull(species)

log_info("Starting genome queries for {length(qc_species)} species")

# Function to query for full genomes
query_full_genome <- function(species_name, query_index, total_queries) {
  log_info("Query {query_index}/{total_queries}: {species_name}")

  result <- tibble(
    species = species_name,
    nuclear_genome = FALSE,
    mitochondrial_genome = FALSE,
    nuclear_count = 0,
    mitochondrial_count = 0,
    nuclear_accessions = NA_character_,
    mitochondrial_accessions = NA_character_,
    error = NA_character_
  )

  tryCatch({
    # Query for nuclear genomes
    query <- paste0(species_name, "[Organism]")
    log_info("  Nuclear query: {nuclear_query}")

    nuclear_results <- rentrez::entrez_search(
      db = "genome",
      term = query,
      retmax = 99000
    )

    if (nuclear_results$count > 0) {
      result$nuclear_genome <- TRUE
      result$nuclear_count <- nuclear_results$count
      result$nuclear_accessions <- paste(nuclear_results$ids, collapse = ",")
      log_success("  Found {nuclear_results$count} nuclear genome(s)")
    } else {
      log_info("  No nuclear genome found")
    }

    # Query for mitochondrial genomes
    mito_query <- paste0(species_name, "[Organism] AND (complete genome[Title] OR complete sequence[Title]) AND (mitochondrion[Title] OR mitochondrial[Title])")
    log_info("  Mitochondrial query: {mito_query}")

    mito_results <- rentrez::entrez_search(
      db = "nucleotide",
      term = mito_query,
      retmax = 99000
    )

    if (mito_results$count > 0) {
      result$mitochondrial_genome <- TRUE
      result$mitochondrial_count <- mito_results$count
      result$mitochondrial_accessions <- paste(mito_results$ids, collapse = ",")
      log_success("  Found {mito_results$count} mitochondrial genome(s)")
    } else {
      log_info("  No mitochondrial genome found")
    }


  }, error = function(e) {
    log_error("  Error querying {species_name}: {e$message}")
    result$error <<- e$message
  })

  return(result)
}

# Query all species
genome_results <- map_df(
  seq_along(qc_species),
  \(i) query_full_genome(qc_species[i], i, length(qc_species)),
  .progress = TRUE
)

# Save results
saveRDS(genome_results, "results/ncbi_genome_results.rds")

# Summary statistics
log_success("Completed!")

