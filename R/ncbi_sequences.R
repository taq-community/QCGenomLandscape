#' Build NCBI search query strings from a BDQC species list and primer map
#'
#' With `batch_size > 1`, species sharing the same `query_marker` are
#' combined into a single OR'd query (`"(Sp1[Organism] OR Sp2[Organism] OR
#' ...) AND marker[Gene]"`), cutting the number of Entrez searches by
#' roughly `batch_size`x. Species are no longer recoverable by parsing the
#' query string in that case -- use the `organism` field NCBI returns in
#' [fetch_ncbi_sequences()]'s results instead (already present in its output).
#'
#' @param species_df Data frame from the BDQC species list csv, already
#'   filtered to `rank == "species"` (must have `species` and `group_en` columns)
#' @param query_primers Data frame from the primers-map csv (must have
#'   `group` and `query_marker` columns)
#' @param batch_size Integer, how many species (sharing the same marker) to
#'   OR together per query, default 1 (one query per species, matching the
#'   original per-species scripts). Keep this well under NCBI's query-length
#'   limits -- 20-50 is a reasonable range in practice.
#' @return Character vector of unique Entrez query strings
#' @export
build_ncbi_queries <- function(species_df, query_primers, batch_size = 1) {
  joined <- species_df |>
    dplyr::mutate(group = tolower(group_en)) |>
    dplyr::left_join(
      query_primers |> dplyr::mutate(group = tolower(group)),
      by = "group"
    ) |>
    dplyr::filter(!is.na(query_marker)) |>
    dplyr::distinct(species, query_marker)

  if (batch_size <= 1) {
    return(
      paste0(joined$species, "[Organism] AND ", joined$query_marker) |>
        unique()
    )
  }

  joined |>
    dplyr::group_by(query_marker) |>
    dplyr::group_map(function(rows, key) {
      purrr::map_chr(chunk_species(rows$species, batch_size), function(sp) {
        paste0(build_organism_clause(sp), " AND ", key$query_marker)
      })
    }) |>
    unlist() |>
    unique() |>
    as.character()
}

#' Fetch NCBI nucleotide records for a set of species/marker queries
#'
#' The batching/error-handling loop behind the NCBI query step. `search_fn`/
#' `summary_fn` are injectable so the batching, high-ID-count flagging, and
#' error-accumulation logic can be unit tested without live network access or
#' an `NCBI_API_KEY`.
#'
#' Queries are no longer restricted server-side to voucher-backed records
#' (the old `AND voucher[Title]` query suffix); instead every record is kept,
#' and `results$is_voucher` flags whether `"voucher"` appears in its title
#' (case-insensitive) -- the same signal the title-filtered query used to
#' rely on, just applied client-side after a single fetch instead of run
#' twice (once filtered, once not).
#'
#' @param queries Character vector of Entrez query strings, e.g. from
#'   [build_ncbi_queries()]
#' @param retmax Integer, per-search max IDs to retrieve, default 5000
#' @param batch_size Integer, `entrez_summary` batch size, default 200
#' @param high_id_threshold Integer, queries returning more IDs than this are
#'   flagged in `high_id_queries`, default 500
#' @param search_fn Function with signature `(db, term, retmax)`, default
#'   [rentrez::entrez_search()]
#' @param summary_fn Function with signature `(db, id)`, default
#'   [rentrez::entrez_summary()]
#' @param progress Logical, show a progress bar, default `TRUE`
#' @return A list with three elements: `results` (tibble of parsed sequence
#'   summaries, with an `is_voucher` logical column), `deficient_queries`
#'   (list of queries that errored), and `high_id_queries` (list of queries
#'   whose ID count exceeded `high_id_threshold`)
#' @importFrom rlang %||%
#' @export
fetch_ncbi_sequences <- function(queries,
                                  retmax = 5000,
                                  batch_size = 200,
                                  high_id_threshold = 500,
                                  search_fn = rentrez::entrez_search,
                                  summary_fn = rentrez::entrez_summary,
                                  progress = TRUE) {
  deficient_queries <- list()
  high_id_queries <- list()

  results <- purrr::map_df(seq_along(queries), \(i) {
    q <- queries[i]
    logger::log_info("Query {i}/{length(queries)}: {q}")

    tryCatch(
      {
        # A single entrez_search(retmax = retmax) already returns both
        # `count` and `ids` in one response -- no need for a separate
        # retmax = 0 call just to read `count` first.
        id_result <- tryCatch(
          search_fn(db = "nucleotide", term = q, retmax = retmax),
          error = function(e) {
            logger::log_error("Error searching for query {i}: {e$message}")
            deficient_queries[[length(deficient_queries) + 1]] <<- list(
              query_index = i,
              query = q,
              error_type = "entrez_search_ids",
              error_message = e$message,
              timestamp = Sys.time()
            )
            NULL
          }
        )

        if (is.null(id_result)) {
          return(tibble::tibble())
        }

        if (id_result$count == 0) {
          logger::log_warn("No results found for query")
          return(tibble::tibble())
        }

        logger::log_info("Found {id_result$count} total sequences, retrieving all IDs...")

        retrieved_ids <- length(id_result$ids)
        logger::log_success("Retrieved {retrieved_ids} IDs")

        if (retrieved_ids > high_id_threshold) {
          high_id_queries[[length(high_id_queries) + 1]] <<- list(
            query_index = i,
            query = q,
            id_count = retrieved_ids,
            timestamp = Sys.time()
          )
          logger::log_warn("High ID count ({retrieved_ids} > {high_id_threshold}) - query stored")
        }

        all_summaries <- list()

        for (batch_start in seq(1, length(id_result$ids), by = batch_size)) {
          batch_end <- min(batch_start + batch_size - 1, length(id_result$ids))
          batch_ids <- id_result$ids[batch_start:batch_end]

          logger::log_info("Fetching summaries {batch_start}-{batch_end} of {length(id_result$ids)}")

          batch_summary <- tryCatch(
            summary_fn(db = "nucleotide", id = batch_ids),
            error = function(e) {
              logger::log_error("Error fetching summaries for batch {batch_start}-{batch_end}: {e$message}")
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

          if (!is.null(batch_summary)) {
            if (length(batch_ids) == 1) {
              all_summaries[[length(all_summaries) + 1]] <- batch_summary
            } else {
              all_summaries <- c(all_summaries, batch_summary)
            }
          }
        }

        purrr::map_df(all_summaries, function(x) {
          if (is.list(x)) {
            subtypes <- tryCatch(strsplit(x$subtype, "\\|")[[1]], error = function(e) character(0))
            subnames <- tryCatch(strsplit(x$subname, "\\|")[[1]], error = function(e) character(0))

            subtype_values <- if (length(subtypes) > 0 && length(subnames) > 0) {
              stats::setNames(as.list(subnames), subtypes)
            } else {
              list()
            }

            tibble::tibble(
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
        }) |> dplyr::mutate(query = q)
      },
      error = function(e) {
        logger::log_error("Error processing query {i}: {e$message}")
        deficient_queries[[length(deficient_queries) + 1]] <<- list(
          query_index = i,
          query = q,
          error_type = "entrez_search",
          error_message = e$message,
          timestamp = Sys.time()
        )
        tibble::tibble()
      }
    )
  }, .progress = progress) |>
    dplyr::filter(!dplyr::if_all(dplyr::everything(), is.na))

  if (nrow(results) > 0) {
    results <- results |>
      dplyr::mutate(is_voucher = dplyr::coalesce(
        grepl("voucher", title, ignore.case = TRUE),
        FALSE
      ))
  }

  list(
    results = results,
    deficient_queries = deficient_queries,
    high_id_queries = high_id_queries
  )
}
