#' Build NCBI search query strings from a BDQC species list and primer map
#'
#' @param species_df Data frame from the BDQC species list csv, already
#'   filtered to `rank == "species"` (must have `species` and `group_en` columns)
#' @param query_primers Data frame from the primers-map csv (must have
#'   `group` and `query_marker` columns)
#' @param voucher Logical; if `TRUE` (default), appends `AND voucher[Title]`
#'   to each query so results are restricted to specimen-voucher-backed records
#' @return Character vector of unique, non-`NA` Entrez query strings
#' @export
build_ncbi_queries <- function(species_df, query_primers, voucher = TRUE) {
  suffix <- if (voucher) " AND voucher[Title]" else ""

  joined <- species_df |>
    dplyr::mutate(group = tolower(group_en)) |>
    dplyr::left_join(
      query_primers |> dplyr::mutate(group = tolower(group)),
      by = "group"
    ) |>
    dplyr::mutate(
      query = ifelse(
        !is.na(query_marker),
        glue::glue("{species}[Organism] AND {query_marker}{suffix}"),
        NA
      )
    )

  joined |>
    dplyr::pull(query) |>
    stats::na.omit() |>
    unique() |>
    as.character()
}

#' Fetch NCBI nucleotide records for a set of species/marker queries
#'
#' The batching/error-handling loop shared by the voucher and non-voucher
#' NCBI query scripts. `search_fn`/`summary_fn` are injectable so the
#' batching, high-ID-count flagging, and error-accumulation logic can be unit
#' tested without live network access or an `NCBI_API_KEY`.
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
#'   summaries), `deficient_queries` (list of queries that errored), and
#'   `high_id_queries` (list of queries whose ID count exceeded `high_id_threshold`)
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
        count_result <- tryCatch(
          search_fn(db = "nucleotide", term = q, retmax = 0),
          error = function(e) {
            logger::log_error("Error in initial search for query {i}: {e$message}")
            deficient_queries[[length(deficient_queries) + 1]] <<- list(
              query_index = i,
              query = q,
              error_type = "entrez_search_count",
              error_message = e$message,
              timestamp = Sys.time()
            )
            NULL
          }
        )

        if (is.null(count_result)) {
          return(tibble::tibble())
        }

        if (count_result$count == 0) {
          logger::log_warn("No results found for query")
          return(tibble::tibble())
        }

        logger::log_info("Found {count_result$count} total sequences, retrieving all IDs...")

        id_result <- tryCatch(
          search_fn(db = "nucleotide", term = q, retmax = retmax),
          error = function(e) {
            logger::log_error("Error retrieving IDs for query {i}: {e$message}")
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

  list(
    results = results,
    deficient_queries = deficient_queries,
    high_id_queries = high_id_queries
  )
}
