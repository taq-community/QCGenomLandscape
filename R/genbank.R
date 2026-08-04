#' Fetch gene annotations for a set of accessions via GenBank XML
#'
#' @param accessions Character vector of NCBI accession numbers
#' @param batch_size Integer, number of accessions fetched per request, default 200
#' @param fetch_fn Function with signature `(db, id, rettype, retmode)`,
#'   default [rentrez::entrez_fetch()]; injectable for testing without network access
#' @param progress Logical, show a progress bar, default `TRUE`
#' @return Tibble with columns `accession`, `gene`, `location`
#' @export
fetch_gene_annotations <- function(accessions, batch_size = 200,
                                    fetch_fn = rentrez::entrez_fetch,
                                    progress = TRUE) {
  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))

  purrr::map_df(seq_along(batches), function(i) {
    batch <- batches[[i]]
    tryCatch(
      {
        gb_raw <- fetch_fn(db = "nucleotide", id = batch, rettype = "gb", retmode = "xml")
        gb <- xml2::read_xml(gb_raw, options = "HUGE")
        seq_nodes <- xml2::xml_find_all(gb, ".//GBSeq")
        purrr::map_df(seq_nodes, function(seq_node) {
          acc <- xml2::xml_text(xml2::xml_find_first(seq_node, ".//GBSeq_accession-version"))
          gene_nodes <- xml2::xml_find_all(seq_node, ".//GBFeature[GBFeature_key='gene']")
          purrr::map_df(gene_nodes, function(node) {
            gene_name <- xml2::xml_text(
              xml2::xml_find_first(node, ".//GBQualifier[GBQualifier_name='gene']/GBQualifier_value")
            )
            location <- xml2::xml_text(xml2::xml_find_first(node, ".//GBFeature_location"))
            tibble::tibble(accession = acc, gene = gene_name, location = location)
          })
        })
      },
      error = function(e) {
        warning(sprintf("Batch %d/%d failed: %s", i, length(batches), e$message))
        tibble::tibble(accession = character(), gene = character(), location = character())
      }
    )
  }, .progress = progress)
}

#' Fetch raw GenBank flat-file text for a set of accessions
#'
#' @param accessions Character vector of NCBI accession numbers
#' @param batch_size Integer, number of accessions fetched per request, default 50
#' @param fetch_fn Function with signature `(db, id, rettype, retmode)`,
#'   default [rentrez::entrez_fetch()]; injectable for testing without network access
#' @param progress Logical, show a progress bar, default `TRUE`
#' @return A list of character scalars, one raw GenBank flat-file blob per batch
#' @export
fetch_gb_records <- function(accessions, batch_size = 50,
                              fetch_fn = rentrez::entrez_fetch,
                              progress = TRUE) {
  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))

  purrr::map(batches, \(batch) {
    fetch_fn(db = "nucleotide", id = batch, rettype = "gb", retmode = "text")
  }, .progress = progress)
}

#' Parse accession/gene/sequence out of GenBank flat-file text
#'
#' @param gb_text Character scalar containing one or more `"//"`-delimited
#'   GenBank flat-file records
#' @return Tibble with columns `accession`, `gene` (`;`-separated if a record
#'   has multiple `/gene=` tags), `sequence` (lowercase-normalized `acgt` only,
#'   `NA` if the record has no `ORIGIN` block)
#' @export
parse_gb_records <- function(gb_text) {
  # Split into individual records at the "//" delimiter
  records <- strsplit(gb_text, "(?m)^//\\s*$", perl = TRUE)[[1]]
  records <- records[nzchar(trimws(records))]

  purrr::map_dfr(records, \(rec) {
    accession <- regmatches(rec, regexpr("(?m)^ACCESSION\\s+(\\S+)", rec, perl = TRUE))
    accession <- if (length(accession)) sub("ACCESSION\\s+", "", accession) else NA_character_

    # Extract gene names from /gene="..." tags
    gene_matches <- regmatches(rec, gregexpr('/gene="([^"]+)"', rec, perl = TRUE))[[1]]
    genes <- unique(sub('/gene="([^"]+)"', "\\1", gene_matches))
    gene <- if (length(genes)) paste(genes, collapse = ";") else NA_character_

    # Extract the ORIGIN sequence block, dropping the "ORIGIN" header line
    # itself before filtering to acgt characters -- otherwise the "G" in
    # "ORIGIN" leaks into the sequence as a spurious leading base
    origin <- regmatches(rec, regexpr("(?ms)^ORIGIN.*", rec, perl = TRUE))
    sequence <- if (length(origin) && nzchar(origin)) {
      origin_body <- sub("^ORIGIN[^\n]*\n", "", origin, perl = TRUE)
      gsub("[^acgtACGT]", "", origin_body)
    } else {
      NA_character_
    }

    tibble::tibble(accession = accession, gene = gene, sequence = sequence)
  })
}
