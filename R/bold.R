#' Fetch BOLD Systems records for a query, optionally saving to disk
#'
#' @param query Character, BOLD query string, default `"geo:province/state:Quebec"`
#' @param extent Character, default `"full"`
#' @param out_path Character or `NULL`; if not `NULL`, writes the raw TSV there
#' @param request_fn Function with signature `(url)` returning an `httr2`
#'   request, default [httr2::request()]; injectable for testing without
#'   network access
#' @return Invisibly, the raw TSV text (character scalar)
#' @export
fetch_bold_sequences <- function(query = "geo:province/state:Quebec",
                                  extent = "full",
                                  out_path = NULL,
                                  request_fn = httr2::request) {
  logger::log_info("BOLD query: {query}")

  query_request <- request_fn("https://portal.boldsystems.org/api/query") |>
    httr2::req_url_query(query = query, extent = extent) |>
    httr2::req_perform()

  query_id <- httr2::resp_body_json(query_request)$query_id
  logger::log_info("BOLD query_id {query_id} -- downloading TSV...")

  download_request <- request_fn(
    paste0("https://portal.boldsystems.org/api/documents/", query_id, "/download")
  ) |>
    httr2::req_url_query(format = "tsv") |>
    httr2::req_perform()

  bold_data <- httr2::resp_body_string(download_request)
  logger::log_success("BOLD download complete ({format(nchar(bold_data), big.mark = ',')} characters)")

  if (!is.null(out_path)) {
    writeLines(bold_data, out_path)
    logger::log_info("BOLD data written to {out_path}")
  }

  invisible(bold_data)
}
