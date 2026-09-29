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
    if (nrow(joined) == 0) {
      return(character(0))
    }
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

#' Retry an Entrez call across transient failures
#'
#' Entrez fails transiently often enough that a single attempt is not a
#' reliable read: DNS resolution times out, the load balancer returns 502 on
#' large batches, and rentrez occasionally throws `"subscript out of bounds"`
#' parsing a truncated response. Every one of those is retryable, and every
#' one of them silently cost this pipeline a whole batch of species before
#' this helper existed -- the 2026-08-19 run lost all 24 Quebec amphibians to
#' a single 10-second DNS timeout on query 1 of 954.
#'
#' @param expr Quoted expression performing the Entrez call
#' @param max_attempts Integer, total attempts before giving up, default 4
#' @param sleep_fn Function with signature `(seconds)`, default [Sys.sleep()];
#'   injectable so retry backoff doesn't slow down tests
#' @param label Character, short description used in log messages
#' @return The value of `expr`, or a `condition` object if every attempt failed
#' @keywords internal
with_entrez_retry <- function(expr, max_attempts = 4, sleep_fn = Sys.sleep,
                               label = "entrez call") {
  expr <- substitute(expr)
  # Capture the caller's frame ONCE, here. Calling parent.frame() inside the
  # tryCatch() below resolves against tryCatch's own stack instead of the
  # caller's, so the expression would be evaluated in the wrong environment.
  caller <- parent.frame()
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(
      eval(expr, envir = caller),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(result)) {
      if (attempt > 1) {
        logger::log_success("{label} succeeded on attempt {attempt}/{max_attempts}")
      }
      return(result)
    }
    logger::log_warn(
      "{label} attempt {attempt}/{max_attempts} failed: {last_error$message}"
    )
    if (attempt < max_attempts) sleep_fn(2^attempt)
  }
  last_error
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
#' @param max_attempts Integer, attempts per Entrez call before the query is
#'   recorded as deficient, default 4 (exponential backoff: 2s, 4s, 8s).
#'   Transient timeouts and 502s are the norm at this query volume, so a
#'   single attempt loses whole batches of species -- see [with_entrez_retry()].
#' @param sleep_fn Function with signature `(seconds)`, default [Sys.sleep()];
#'   injectable so retry backoff doesn't slow down tests
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
                                  max_attempts = 4,
                                  sleep_fn = Sys.sleep,
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
        id_result <- with_entrez_retry(
          search_fn(db = "nucleotide", term = q, retmax = retmax),
          max_attempts = max_attempts,
          sleep_fn = sleep_fn,
          label = paste0("search query ", i)
        )

        if (inherits(id_result, "condition")) {
          logger::log_error("Error searching for query {i}: {id_result$message}")
          deficient_queries[[length(deficient_queries) + 1]] <<- list(
            query_index = i,
            query = q,
            error_type = "entrez_search_ids",
            error_message = id_result$message,
            timestamp = Sys.time()
          )
          return(tibble::tibble())
        }

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

          batch_summary <- with_entrez_retry(
            summary_fn(db = "nucleotide", id = batch_ids),
            max_attempts = max_attempts,
            sleep_fn = sleep_fn,
            label = paste0("summaries ", batch_start, "-", batch_end, " of query ", i)
          )

          if (inherits(batch_summary, "condition")) {
            logger::log_error("Error fetching summaries for batch {batch_start}-{batch_end}: {batch_summary$message}")
            deficient_queries[[length(deficient_queries) + 1]] <<- list(
              query_index = i,
              query = q,
              error_type = "entrez_summary",
              error_message = batch_summary$message,
              batch_range = paste0(batch_start, "-", batch_end),
              batch_ids = batch_ids,
              timestamp = Sys.time()
            )
            batch_summary <- NULL
          }

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

#' Summarise the species lost to failed Entrez queries
#'
#' `fetch_ncbi_sequences()` records every query it could not complete in
#' `deficient_queries` and carries on, so a failed batch does not abort a
#' multi-hour run. Nothing downstream consumed that record, which is how a
#' single timed-out batch became a published "Amphibians: 0% coverage" panel.
#' This turns the failure log back into the list of species it cost, so the
#' pipeline can surface -- or refuse to pass -- an incomplete corpus.
#'
#' @param deficient_queries List as returned in
#'   `fetch_ncbi_sequences()$deficient_queries`
#' @param species Character vector of species names to look for in the failed
#'   query strings, typically the BDQC species list
#' @return Tibble with one row per (query_index, species) lost, columns
#'   `query_index`, `error_type`, `error_message`, `species`
#' @export
deficient_query_species <- function(deficient_queries, species) {
  empty <- tibble::tibble(
    query_index = integer(), error_type = character(),
    error_message = character(), species = character()
  )
  if (!length(deficient_queries)) {
    return(empty)
  }
  purrr::map_df(deficient_queries, function(d) {
    hit <- species[vapply(species, grepl, logical(1), x = d$query, fixed = TRUE)]
    if (!length(hit)) {
      return(empty)
    }
    tibble::tibble(
      query_index = d$query_index,
      error_type = d$error_type,
      error_message = d$error_message,
      species = hit
    )
  })
}

#' Stop the pipeline when Entrez queries were lost
#'
#' Called as its own `targets` target so that an incomplete NCBI corpus is a
#' build failure with the affected species named, rather than a silent hole
#' that only shows up as an implausible zero in a figure. Set
#' `allow_deficient = TRUE` (or `QCGENOM_ALLOW_DEFICIENT=true`) to downgrade
#' it to a warning when a partial corpus is genuinely acceptable.
#'
#' @param deficient_queries List as returned in
#'   `fetch_ncbi_sequences()$deficient_queries`
#' @param species Character vector of species names, typically the BDQC list
#' @param allow_deficient Logical, warn instead of erroring
#' @return The tibble from [deficient_query_species()], invisibly
#' @export
assert_no_deficient_queries <- function(deficient_queries, species,
                                        allow_deficient = identical(
                                          tolower(Sys.getenv("QCGENOM_ALLOW_DEFICIENT")), "true"
                                        )) {
  lost <- deficient_query_species(deficient_queries, species)
  if (nrow(lost) == 0) {
    return(invisible(lost))
  }
  msg <- sprintf(
    paste0(
      "%d Entrez quer%s failed after retries, costing %d species.\n",
      "Affected species are absent from the corpus and will read as zero coverage.\n",
      "Errors: %s\nFirst species lost: %s"
    ),
    length(unique(lost$query_index)),
    if (length(unique(lost$query_index)) == 1) "y" else "ies",
    length(unique(lost$species)),
    paste(unique(lost$error_type), collapse = ", "),
    paste(utils::head(unique(lost$species), 8), collapse = ", ")
  )
  if (allow_deficient) {
    warning(msg, call. = FALSE)
  } else {
    stop(msg, call. = FALSE)
  }
  invisible(lost)
}

#' Drop subspecies/infraspecific records from NCBI results
#'
#' NCBI's `[Organism]` search matches a species name *and* everything below
#' it taxonomically (subspecies, breeds, domestic forms), e.g. querying
#' `"Canis lupus[Organism]"` also returns `"Canis lupus familiaris"` records.
#' Left unfiltered, these can dominate a species' sequence count (domestic
#' dog barcodes outnumber wild wolf ones ~3:1 in the raw NCBI pull) and bias
#' any coverage/quality analysis keyed on the BDQC species list. Keeping only
#' rows whose `organism` exactly matches a name in `species` removes them.
#'
#' @param results Tibble with an `organism` column, e.g. `fetch_ncbi_sequences()$results`
#' @param species Character vector of exact species names to keep (typically
#'   the BDQC species list)
#' @return `results`, filtered to rows where `organism %in% species`
#' @export
filter_named_species <- function(results, species) {
  n_before <- nrow(results)
  filtered <- dplyr::filter(results, organism %in% species)
  logger::log_info(
    "Subspecies filter: kept {nrow(filtered)}/{n_before} records ({n_before - nrow(filtered)} infraspecific records dropped)"
  )
  filtered
}
