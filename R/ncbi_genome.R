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

#' Resolve nuclear genome availability for one batch of species
#'
#' Nuclear genome assemblies are rare for non-model species -- an OR'd
#' existence check (`retmax = 0`) across the whole batch almost always comes
#' back empty, so most batches skip straight to "no genome" from a single
#' call. Only a batch with an actual hit pays for per-species resolution.
#' `genome`'s `entrez_summary()` schema isn't reliably documented, so unlike
#' `resolve_mitochondrial_genomes()` this can't attribute a batched hit back
#' to species that way -- it re-queries individually instead.
#'
#' @param batch_species Character vector of species in this batch
#' @param organism_clause The batch's OR'd `[Organism]` clause, from
#'   `build_organism_clause()`
#' @param search_fn Function with signature `(db, term, retmax)`
#' @return Tibble with columns `species`, `nuclear_genome`, `nuclear_count`,
#'   `nuclear_accessions`
#' @noRd
resolve_nuclear_genomes <- function(batch_species, organism_clause, search_fn) {
  hit <- tryCatch(
    search_fn(db = "genome", term = organism_clause, retmax = 0)$count > 0,
    error = function(e) TRUE # existence check itself failed -- fall back to per-species, don't assume "no genome"
  )

  if (!hit) {
    return(tibble::tibble(
      species = batch_species, nuclear_genome = FALSE,
      nuclear_count = 0, nuclear_accessions = NA_character_
    ))
  }

  purrr::map_df(batch_species, function(sp) {
    res <- tryCatch(
      search_fn(db = "genome", term = build_organism_clause(sp), retmax = 99000),
      error = function(e) list(count = 0, ids = character(0))
    )
    tibble::tibble(
      species = sp,
      nuclear_genome = res$count > 0,
      nuclear_count = res$count,
      nuclear_accessions = if (res$count > 0) paste(res$ids, collapse = ",") else NA_character_
    )
  })
}

#' Resolve mitochondrial genome availability for one batch of species
#'
#' Unlike nuclear assemblies, "complete mitochondrial genome/sequence"-titled
#' records are common enough that an existence pre-check rarely saves
#' anything -- it would just add an extra round trip before doing the same
#' work anyway. So this goes straight to one batched search across the whole
#' batch, then attributes each hit back to its species via
#' `fetch_summaries_batched()`'s `organism` field -- the nucleotide
#' database's summary schema reliably includes it (same field
#' [fetch_ncbi_sequences()] already relies on).
#'
#' @inheritParams resolve_nuclear_genomes
#' @param summary_fn Function with signature `(db, id)`
#' @return Tibble with columns `species`, `mitochondrial_genome`,
#'   `mitochondrial_count`, `mitochondrial_accessions`
#' @noRd
resolve_mitochondrial_genomes <- function(batch_species, organism_clause, search_fn, summary_fn) {
  empty <- tibble::tibble(
    species = batch_species, mitochondrial_genome = FALSE,
    mitochondrial_count = 0, mitochondrial_accessions = NA_character_
  )

  mito_query <- paste(
    organism_clause,
    "AND (complete genome[Title] OR complete sequence[Title]) AND (mitochondrion[Title] OR mitochondrial[Title])"
  )
  result <- tryCatch(search_fn(db = "nucleotide", term = mito_query, retmax = 99000), error = function(e) NULL)

  if (is.null(result) || result$count == 0 || length(result$ids) == 0) {
    return(empty)
  }

  summaries <- fetch_summaries_batched(result$ids, summary_fn = summary_fn, db = "nucleotide")

  hits <- purrr::map_df(summaries, function(x) {
    tibble::tibble(id = x$uid %||% NA_character_, species = x$organism %||% NA_character_)
  })
  hits <- hits[hits$species %in% batch_species, , drop = FALSE]

  if (nrow(hits) == 0) {
    return(empty)
  }

  resolved <- hits |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      mitochondrial_genome = TRUE,
      mitochondrial_count = dplyr::n(),
      mitochondrial_accessions = paste(id, collapse = ","),
      .groups = "drop"
    )

  empty |>
    dplyr::select(species) |>
    dplyr::left_join(resolved, by = "species") |>
    dplyr::mutate(
      mitochondrial_genome = !is.na(mitochondrial_genome) & mitochondrial_genome,
      mitochondrial_count = dplyr::coalesce(mitochondrial_count, 0)
    )
}

#' Query full-genome availability for a vector of species, batched
#'
#' Splits `species` into batches of `batch_size` and resolves nuclear and
#' mitochondrial genome availability independently per batch -- see
#' `resolve_nuclear_genomes()` and `resolve_mitochondrial_genomes()` for why
#' they use different strategies (nuclear hits are rare enough for a fast
#' existence-check path to pay off; mitochondrial hits are not, so that
#' dimension goes straight to a batched search + summary attribution).
#'
#' @param species Character vector of species names (deduplicated internally)
#' @param batch_size Integer, species OR'd together per query, default 25.
#'   Keep well under NCBI's query-length limits.
#' @param search_fn Function with signature `(db, term, retmax)`, default
#'   [rentrez::entrez_search()]; injectable for testing without network access
#' @param summary_fn Function with signature `(db, id)`, default
#'   [rentrez::entrez_summary()]
#' @param progress Logical, show a progress bar, default `TRUE`
#' @return Tibble, one row per species: `species`, `nuclear_genome`,
#'   `mitochondrial_genome`, `nuclear_count`, `mitochondrial_count`,
#'   `nuclear_accessions`, `mitochondrial_accessions`, `error`
#' @export
fetch_ncbi_genomes <- function(species, batch_size = 25,
                                search_fn = rentrez::entrez_search,
                                summary_fn = rentrez::entrez_summary,
                                progress = TRUE) {
  species <- unique(species)
  batches <- chunk_species(species, batch_size)

  purrr::map_df(seq_along(batches), function(b) {
    batch_species <- batches[[b]]
    logger::log_info("Genome batch {b}/{length(batches)}: {length(batch_species)} species")
    organism_clause <- build_organism_clause(batch_species)

    nuclear_df <- resolve_nuclear_genomes(batch_species, organism_clause, search_fn)
    mito_df <- resolve_mitochondrial_genomes(batch_species, organism_clause, search_fn, summary_fn)

    dplyr::full_join(nuclear_df, mito_df, by = "species") |>
      dplyr::mutate(error = NA_character_)
  }, .progress = progress)
}
