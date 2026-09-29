#' Load and normalize a conservation-risk-status table
#'
#' Unifies the CA/COSEPAC and QC/LEMV risk-status loading blocks previously
#' duplicated between `create_dataframe.R` and `taxon_representation.R`.
#'
#' @param path Character, CSV path (`data/CA_especes_en_peril.csv` for `"CA"`,
#'   `data/QC_especes_en_peril.csv` for `"QC"`)
#' @param jurisdiction `"CA"` or `"QC"`
#' @param translate Logical; if `TRUE`, recodes French status labels to
#'   English (as `taxon_representation.R` did); default `FALSE` keeps the
#'   original French labels (as `create_dataframe.R` did)
#' @return Tibble with columns `species`, `status`, `jurisdiction`
#' @export
load_risk_status <- function(path, jurisdiction = c("CA", "QC"), translate = FALSE) {
  jurisdiction <- match.arg(jurisdiction)
  raw <- utils::read.csv(path, fileEncoding = "UTF-8-BOM")

  risk <- if (jurisdiction == "CA") {
    raw |>
      dplyr::mutate(
        # regmatches() silently returns a SHORTER vector when a row does not
        # match, which makes mutate() fail with a size error rather than
        # flagging the offending row. The COSEWIC list contains at least one
        # non-binomial entry -- "Ambystoma (2) laterale - jeffersonianum",
        # the unisexual Ambystoma complex -- so extract length-safely and let
        # the row carry NA, which the binomial filter below then drops.
        species = extract_binomial(`Nom.scientifique`),
        status = `Statut.selon.le.COSEPAC`
      ) |>
      dplyr::filter(status != "" & status != "Non active")
  } else {
    raw |>
      dplyr::mutate(
        species = paste(GENRE, ESPECE),
        status = STATUT_LEMV
      ) |>
      dplyr::filter(status != "Retir\u00e9e" & status != "Non suivie")
  }

  risk <- risk |>
    dplyr::mutate(jurisdiction = jurisdiction) |>
    dplyr::filter(!is.na(species)) |>
    dplyr::select(species, status, jurisdiction) |>
    dplyr::distinct()

  if (translate) {
    risk <- risk |>
      dplyr::mutate(status = dplyr::recode(status,
        "Disparue" = "Extinct",
        "Disparue du pays" = "Extirpated",
        "Donn\u00e9es insuffisantes" = "Data deficient",
        "En voie de disparition" = "Endangered",
        "Menac\u00e9e" = "Threatened",
        "Non en p\u00e9ril" = "Not at risk",
        "Pr\u00e9occupante" = "Special concern",
        "Candidate" = "Candidate",
        "Susceptible" = "Likely to be designated",
        "Vuln\u00e9rable" = "Vulnerable"
      ))
  }

  risk
}

#' Extract a `Genus species` binomial, preserving vector length
#'
#' [regmatches()] drops non-matching elements, so using it inside
#' [dplyr::mutate()] turns one malformed name into a whole-pipeline size
#' error. This returns `NA` for rows that carry no binomial instead.
#'
#' @param x Character vector of scientific names
#' @return Character vector the same length as `x`
#' @examples
#' extract_binomial(c("Ambystoma maculatum", "Ambystoma (2) laterale - jeffersonianum"))
#' @export
extract_binomial <- function(x) {
  m <- regexpr("^\\w+\\s+\\w+", x)
  out <- rep(NA_character_, length(x))
  out[m != -1] <- regmatches(x, m)
  out
}
