#' Query NCBI genome/nucleotide databases for full-genome availability
#'
#' @param species_name Character scalar, species name
#' @param query_index Integer, for logging only, default 1
#' @param total_queries Integer, for logging only, default 1
#' @param search_fn Function with signature `(db, term, retmax)`, default
#'   [rentrez::entrez_search()]; injectable for testing without network access
#' @return Tibble with one row: `species`, `nuclear_genome`,
#'   `mitochondrial_genome`, `nuclear_count`, `mitochondrial_count`,
#'   `nuclear_accessions`, `mitochondrial_accessions`, `error`
#' @export
query_full_genome <- function(species_name, query_index = 1L, total_queries = 1L,
                               search_fn = rentrez::entrez_search) {
  logger::log_info("Query {query_index}/{total_queries}: {species_name}")

  result <- tibble::tibble(
    species = species_name,
    nuclear_genome = FALSE,
    mitochondrial_genome = FALSE,
    nuclear_count = 0,
    mitochondrial_count = 0,
    nuclear_accessions = NA_character_,
    mitochondrial_accessions = NA_character_,
    error = NA_character_
  )

  tryCatch(
    {
      query <- paste0(species_name, "[Organism]")
      logger::log_info("  Nuclear query: {query}")

      nuclear_results <- search_fn(db = "genome", term = query, retmax = 99000)

      if (nuclear_results$count > 0) {
        result$nuclear_genome <- TRUE
        result$nuclear_count <- nuclear_results$count
        result$nuclear_accessions <- paste(nuclear_results$ids, collapse = ",")
        logger::log_success("  Found {nuclear_results$count} nuclear genome(s)")
      } else {
        logger::log_info("  No nuclear genome found")
      }

      mito_query <- paste0(
        species_name,
        "[Organism] AND (complete genome[Title] OR complete sequence[Title]) AND (mitochondrion[Title] OR mitochondrial[Title])"
      )
      logger::log_info("  Mitochondrial query: {mito_query}")

      mito_results <- search_fn(db = "nucleotide", term = mito_query, retmax = 99000)

      if (mito_results$count > 0) {
        result$mitochondrial_genome <- TRUE
        result$mitochondrial_count <- mito_results$count
        result$mitochondrial_accessions <- paste(mito_results$ids, collapse = ",")
        logger::log_success("  Found {mito_results$count} mitochondrial genome(s)")
      } else {
        logger::log_info("  No mitochondrial genome found")
      }
    },
    error = function(e) {
      logger::log_error("  Error querying {species_name}: {e$message}")
      result$error <<- e$message
    }
  )

  result
}

#' Query full-genome availability for a vector of species
#'
#' @param species Character vector of species names
#' @param progress Logical, show a progress bar, default `TRUE`
#' @param ... Passed through to [query_full_genome()] (e.g. `search_fn` for testing)
#' @return Tibble, one row per species (see [query_full_genome()])
#' @export
fetch_ncbi_genomes <- function(species, progress = TRUE, ...) {
  purrr::map_df(
    seq_along(species),
    \(i) query_full_genome(species[i], i, length(species), ...),
    .progress = progress
  )
}
