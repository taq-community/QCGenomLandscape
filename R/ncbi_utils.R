#' Build an OR'd Entrez `[Organism]` clause from a vector of species names
#'
#' Shared by [build_ncbi_queries()]'s batched branch and [fetch_ncbi_genomes()]
#' -- both need to combine multiple species into one Entrez search term.
#'
#' @param species Character vector of species names
#' @return Character scalar, e.g. `"(Sp1[Organism] OR Sp2[Organism])"`
#' @noRd
build_organism_clause <- function(species) {
  paste0("(", paste(paste0(species, "[Organism]"), collapse = " OR "), ")")
}

#' Split a vector into `batch_size`-sized chunks, in order
#'
#' Shared by [build_ncbi_queries()]'s batched branch and [fetch_ncbi_genomes()].
#'
#' @param x A vector
#' @param batch_size Integer, max chunk size
#' @return A list of vectors, each of length at most `batch_size`
#' @noRd
chunk_species <- function(x, batch_size) {
  split(x, ceiling(seq_along(x) / batch_size))
}

#' Fetch `entrez_summary()` records in ID batches, flattened to one list
#'
#' `entrez_summary()` can hit an HTTP 414 if given too many IDs at once, and
#' -- confusingly -- returns an unwrapped single record (not a length-1 list
#' containing one record) when given exactly one ID. Both quirks need
#' handling wherever summaries are fetched for a set of IDs; used by the
#' mitochondrial-genome resolution in [fetch_ncbi_genomes()].
#'
#' @param ids Character vector of NCBI UIDs
#' @param summary_fn Function with signature `(db, id)`, default
#'   [rentrez::entrez_summary()]
#' @param db Character, database name, default `"nucleotide"`
#' @param batch_size Integer, IDs per `entrez_summary()` call, default 200
#' @return A list of summary records (each itself a list), one per ID that
#'   didn't error
#' @noRd
fetch_summaries_batched <- function(ids, summary_fn = rentrez::entrez_summary, db = "nucleotide", batch_size = 200) {
  all_summaries <- list()

  for (batch_start in seq(1, length(ids), by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, length(ids))
    batch_ids <- ids[batch_start:batch_end]

    batch_summary <- tryCatch(summary_fn(db = db, id = batch_ids), error = function(e) NULL)

    if (!is.null(batch_summary)) {
      if (length(batch_ids) == 1) {
        all_summaries[[length(all_summaries) + 1]] <- batch_summary
      } else {
        all_summaries <- c(all_summaries, batch_summary)
      }
    }
  }

  all_summaries
}
