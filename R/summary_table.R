#' Build the per-species genomic-data summary table
#'
#' Pure, in-memory core of `create_dataframe.R`: joins NCBI sequence counts,
#' gene-marker counts, BDQC taxonomy, and conservation-risk status into one
#' row per species.
#'
#' @param ncbi_results Tibble as returned by `fetch_ncbi_sequences()$results`
#'   (must have `organism`, `accession` columns)
#' @param genes_df Tibble as returned by [fetch_gene_annotations()] (must
#'   have `accession`, `gene` columns)
#' @param bdqc_taxo Tibble with columns `species`, `vernacular_fr`,
#'   `vernacular_en`, `group_en`
#' @param ca_risk,qc_risk Tibble with columns `species`, `status` (e.g. from
#'   [load_risk_status()])
#' @return Tibble, one row per species
#' @export
build_summary_dataframe <- function(ncbi_results, genes_df, bdqc_taxo, ca_risk, qc_risk) {
  # `organism` is NCBI's own per-record taxon name -- more reliable than
  # parsing it back out of `query`, and required once build_ncbi_queries()
  # batches multiple species into one OR'd query (query no longer starts
  # with a single species name in that case).
  ncbi_results <- ncbi_results |>
    dplyr::mutate(species = organism)

  gene_counts <- genes_df |>
    dplyr::mutate(gene_group = assign_gene_group(gene)) |>
    dplyr::left_join(ncbi_results |> dplyr::select(accession, species), by = "accession") |>
    dplyr::filter(!is.na(species), gene_group != "Other") |>
    dplyr::distinct(species, accession, gene_group) |>
    dplyr::count(species, gene_group, name = "n_seq") |>
    tidyr::pivot_wider(names_from = gene_group, values_from = n_seq, values_fill = 0) |>
    dplyr::distinct()

  total_seq <- ncbi_results |>
    dplyr::count(species, name = "n_total_seq")

  # Column names in c() are quoted strings, not backtick symbols, so \uxxxx
  # escapes work here -- keeps this source file ASCII-only per R CMD check.
  display_names <- c(
    "Nom scientifique" = "species",
    "Nom commun FR" = "vernacular_fr",
    "Nom commun EN" = "vernacular_en",
    "Groupe taxonomique" = "group_en",
    "Statut de l'esp\u00e8ce au Canada" = "statut_canada",
    "Statut de l'esp\u00e8ce au Qu\u00e9bec" = "statut_quebec",
    "S\u00e9quences totales (NCBI)" = "n_total_seq"
  )

  ncbi_results |>
    dplyr::distinct(species) |>
    dplyr::filter(!is.na(species)) |>
    dplyr::left_join(bdqc_taxo, by = "species") |>
    dplyr::left_join(ca_risk |> dplyr::select(species, statut_canada = status), by = "species") |>
    dplyr::left_join(qc_risk |> dplyr::select(species, statut_quebec = status), by = "species") |>
    dplyr::left_join(total_seq, by = "species") |>
    dplyr::left_join(gene_counts, by = "species") |>
    dplyr::select(
      species, vernacular_fr, vernacular_en, group_en,
      statut_canada, statut_quebec, n_total_seq,
      dplyr::any_of(c(
        "COI", "Cytb", "ND1", "ND2", "ND4", "ND5",
        "12S rRNA", "16S rRNA", "Nuclear rRNA / ITS",
        "Photosynthesis-related (rbcL, matK, etc.)"
      ))
    ) |>
    dplyr::rename(!!!display_names)
}

#' Build the per-species genomic-data summary table, reading inputs from disk
#'
#' Thin I/O wrapper around [build_summary_dataframe()] for interactive/ad hoc
#' use outside the `_targets.R` pipeline (which calls [build_summary_dataframe()]
#' directly with in-memory pipeline objects).
#'
#' @param ncbi_path Character, path to `fetch_ncbi_sequences()$results`, saved as RDS
#' @param genes_path Character, path to [fetch_gene_annotations()] output, saved as RDS
#' @param bdqc_path Character, path to the BDQC species list CSV
#' @param ca_risk_path,qc_risk_path Character, paths to the CA/QC risk-status CSVs
#' @return Tibble, one row per species (see [build_summary_dataframe()])
#' @export
build_summary_dataframe_from_files <- function(
    ncbi_path = "results/ncbi_results.rds",
    genes_path = "results/genes_subsamp_50_df.rds",
    bdqc_path = "data/bdqc_list_01122025.csv",
    ca_risk_path = "data/CA_especes_en_peril.csv",
    qc_risk_path = "data/QC_especes_en_peril.csv") {
  ncbi_results <- readRDS(ncbi_path)
  genes_df <- readRDS(genes_path)

  bdqc_taxo <- utils::read.csv(bdqc_path) |>
    dplyr::filter(rank == "species") |>
    dplyr::select(species, vernacular_fr, vernacular_en, group_en) |>
    dplyr::distinct()

  ca_risk <- load_risk_status(ca_risk_path, jurisdiction = "CA")
  qc_risk <- load_risk_status(qc_risk_path, jurisdiction = "QC")

  build_summary_dataframe(ncbi_results, genes_df, bdqc_taxo, ca_risk, qc_risk)
}
