#' @keywords internal
#' @importFrom rlang .data .env
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Column names referenced unquoted inside dplyr/case_when NSE calls --
# not undefined globals, just data-frame columns R CMD check can't see
# statically.
utils::globalVariables(c(
  "ESPECE", "GENRE", "Nom.scientifique", "STATUT_LEMV", "Statut.selon.le.COSEPAC",
  "accession", "gene", "gene_group", "group", "group_en", "groupe", "has_data", "id",
  "mitochondrial_count", "mitochondrial_genome", "n_seq", "n_sp", "n_sp_with_data",
  "n_species", "n_total", "n_total_seq", "organism", "pct", "query", "query_marker",
  "species", "status", "statut_canada", "statut_quebec", "total", "vernacular_en",
  "vernacular_fr"
))
