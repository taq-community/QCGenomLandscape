library(targets)
library(tarchetypes)
library(QCGenomLandscape)

tar_option_set(
  packages = c("QCGenomLandscape", "sf", "Biostrings"),
  format = "rds"
)

# ---- One-time setup: logger + NCBI API key -------------------------------

log_dir <- "logs"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

logger::log_appender(logger::appender_file(
  file.path(log_dir, sprintf("targets_pipeline_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
))
logger::log_threshold(logger::INFO)
rentrez::set_entrez_key(Sys.getenv("NCBI_API_KEY"))

# ---- Dev/test overrides ----------------------------------------------------
# QCGL_RESULTS_DIR: write sink targets somewhere other than results/, so a
# test run can never overwrite a real one.
# QCGL_MAX_SPECIES: cap bdqc_species to the first N rows, for a fast smoke
# test of the whole pipeline. Both unset by default (full run, results/).
results_dir <- Sys.getenv("QCGL_RESULTS_DIR", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

max_species <- suppressWarnings(as.integer(Sys.getenv("QCGL_MAX_SPECIES", NA)))
if (!is.na(max_species)) {
  logger::log_info("QCGL_MAX_SPECIES set -- capping to {max_species} species")
}

# ---- Swift export helpers --------------------------------------------------
# Uploads are skipped (with a warning, not an error) if SWIFT_APP_CRED_ID/
# SWIFT_APP_CRED_SECRET aren't set, so the pipeline still runs end-to-end
# without Swift configured. Application credentials, not username/password,
# because Arbutus accounts with MFA can't authenticate with a plain password
# over the API -- see ?swift_auth for how to generate one.

get_swift_token <- function() {
  required <- c("SWIFT_APP_CRED_ID", "SWIFT_APP_CRED_SECRET")
  if (!all(nzchar(Sys.getenv(required)))) {
    msg <- "SWIFT_APP_CRED_ID/SWIFT_APP_CRED_SECRET not set -- Swift upload targets will be skipped"
    logger::log_warn(msg)
    warning(msg, call. = FALSE)
    return(NULL)
  }
  swift_auth()
}

maybe_upload_to_swift <- function(path, token) {
  if (is.null(token)) {
    logger::log_warn("Skipping Swift upload of {path} (no token)")
    warning("Skipping Swift upload of ", path, " (no token)", call. = FALSE)
    return(invisible(NULL))
  }
  upload_to_swift(path, token = token)
}

list(
  # ---- Inputs ----------------------------------------------------------
  tar_target(bdqc_csv, "data/bdqc_list_01122025.csv", format = "file"),
  tar_target(primers_csv, "data/primers_map_group_bdqc_list_01122025.csv", format = "file"),
  tar_target(ca_risk_csv, "data/CA_especes_en_peril.csv", format = "file"),
  tar_target(qc_risk_csv, "data/QC_especes_en_peril.csv", format = "file"),
  tar_target(canvec_gdb, "data/canvec_1M_CA_Admin.gdb", format = "file"),

  tar_target(bdqc_species, {
    species <- read.csv(bdqc_csv) |> dplyr::filter(rank == "species")
    if (!is.na(max_species)) species <- head(species, max_species)
    species
  }),
  tar_target(query_primers, read.csv2(primers_csv)),
  tar_target(bdqc_taxo, bdqc_species |>
    dplyr::select(species, vernacular_fr, vernacular_en, group_en) |>
    dplyr::distinct()),
  tar_target(bdqc_grouped, bdqc_species |>
    dplyr::mutate(groupe = QCGenomLandscape::classify_taxon_group(group_en, order, class, phylum, kingdom)) |>
    dplyr::select(species, groupe)),

  tar_target(ca_risk, QCGenomLandscape::load_risk_status(ca_risk_csv, jurisdiction = "CA")),
  tar_target(qc_risk, QCGenomLandscape::load_risk_status(qc_risk_csv, jurisdiction = "QC")),
  tar_target(ca_risk_en, QCGenomLandscape::load_risk_status(ca_risk_csv, jurisdiction = "CA", translate = TRUE)),
  tar_target(qc_risk_en, QCGenomLandscape::load_risk_status(qc_risk_csv, jurisdiction = "QC", translate = TRUE)),

  tar_target(qc_boundary, QCGenomLandscape::load_canvec_boundary(canvec_gdb, jurisdiction = 102)),
  tar_target(ca_boundary, QCGenomLandscape::load_canvec_boundary(canvec_gdb, country = 140)),

  # ---- 1. BOLD -----------------------------------------------------------
  tar_target(bold_raw, QCGenomLandscape::fetch_bold_sequences()),
  tar_target(bold_saved, {
    writeLines(bold_raw, file.path(results_dir, "bold_qc_data.tsv"))
    file.path(results_dir, "bold_qc_data.tsv")
  }, format = "file"),

  # ---- 2. NCBI sequences (voucher) ---------------------------------------
  # batch_size batches species sharing a marker into one OR'd Entrez query --
  # cuts total NCBI requests by roughly that factor. Species are recovered
  # afterwards from the `organism` field NCBI returns, not by parsing `query`.
  tar_target(
    ncbi_queries_voucher,
    QCGenomLandscape::build_ncbi_queries(bdqc_species, query_primers, voucher = TRUE, batch_size = 25)
  ),
  tar_target(ncbi_voucher, QCGenomLandscape::fetch_ncbi_sequences(ncbi_queries_voucher)),
  tar_target(ncbi_results_saved, {
    saveRDS(ncbi_voucher$results, file.path(results_dir, "ncbi_results.rds"))
    file.path(results_dir, "ncbi_results.rds")
  }, format = "file"),
  tar_target(deficient_queries_saved, {
    saveRDS(ncbi_voucher$deficient_queries, file.path(results_dir, "deficient_queries.rds"))
    file.path(results_dir, "deficient_queries.rds")
  }, format = "file"),
  tar_target(high_id_queries_saved, {
    saveRDS(ncbi_voucher$high_id_queries, file.path(results_dir, "high_id_queries.rds"))
    file.path(results_dir, "high_id_queries.rds")
  }, format = "file"),

  # ---- 3. NCBI sequences (non-voucher) + geo-filter ----------------------
  tar_target(
    ncbi_queries_non_voucher,
    QCGenomLandscape::build_ncbi_queries(bdqc_species, query_primers, voucher = FALSE, batch_size = 25)
  ),
  tar_target(ncbi_non_voucher, QCGenomLandscape::fetch_ncbi_sequences(ncbi_queries_non_voucher)),
  tar_target(ncbi_non_voucher_results_saved, {
    saveRDS(ncbi_non_voucher$results, file.path(results_dir, "ncbi_non_voucher_results.rds"))
    file.path(results_dir, "ncbi_non_voucher_results.rds")
  }, format = "file"),
  tar_target(deficient_queries_non_voucher_saved, {
    saveRDS(ncbi_non_voucher$deficient_queries, file.path(results_dir, "deficient_queries_non_voucher_results.rds"))
    file.path(results_dir, "deficient_queries_non_voucher_results.rds")
  }, format = "file"),
  tar_target(high_id_queries_non_voucher_saved, {
    saveRDS(ncbi_non_voucher$high_id_queries, file.path(results_dir, "high_id_queries_non_voucher_results.rds"))
    file.path(results_dir, "high_id_queries_non_voucher_results.rds")
  }, format = "file"),
  # Geo-filtered (QC/CA) view of the non-voucher results -- demonstrates
  # load_canvec_boundary()/flag_within_boundary(), not saved to a named
  # results/ file since the original script only used this interactively
  tar_target(ncbi_non_voucher_geo, {
    pts <- ncbi_non_voucher$results |>
      dplyr::mutate(
        parsed_coords = purrr::map(lat_lon, QCGenomLandscape::parse_latlon),
        latitude = purrr::map_dbl(parsed_coords, "lat"),
        longitude = purrr::map_dbl(parsed_coords, "lon")
      ) |>
      dplyr::select(-parsed_coords) |>
      dplyr::filter(!is.na(latitude) & !is.na(longitude)) |>
      sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

    pts <- QCGenomLandscape::flag_within_boundary(pts, qc_boundary, "in_qc")
    QCGenomLandscape::flag_within_boundary(pts, ca_boundary, "in_ca")
  }),

  # ---- 4. NCBI full genomes -----------------------------------------------
  # Batched: a cheap OR'd existence check per batch of 25 species short-
  # circuits the (common) all-zero-hits case, falling back to per-species
  # resolution only for batches with an actual genome/mitogenome hit.
  tar_target(
    ncbi_genomes,
    QCGenomLandscape::fetch_ncbi_genomes(bdqc_species$species, batch_size = 25)
  ),
  tar_target(ncbi_genome_results_saved, {
    saveRDS(ncbi_genomes, file.path(results_dir, "ncbi_genome_results.rds"))
    file.path(results_dir, "ncbi_genome_results.rds")
  }, format = "file"),

  # ---- 5. Gene annotations (subsample) + taxon_representation figures ----
  tar_target(gene_accessions, {
    set.seed(42)
    ncbi_voucher$results |>
      dplyr::mutate(species = organism) |>
      dplyr::group_by(species) |>
      dplyr::slice_sample(n = 5) |>
      dplyr::ungroup() |>
      dplyr::pull(accession)
  }),
  tar_target(gene_annotations, QCGenomLandscape::fetch_gene_annotations(gene_accessions)),
  tar_target(genes_saved, {
    saveRDS(gene_annotations, file.path(results_dir, "genes_subsamp_50_df.rds"))
    file.path(results_dir, "genes_subsamp_50_df.rds")
  }, format = "file"),

  tar_target(genes_grouped, gene_annotations |>
    dplyr::mutate(gene_group = QCGenomLandscape::assign_gene_group(gene)) |>
    dplyr::left_join(
      ncbi_voucher$results |>
        dplyr::mutate(species = organism) |>
        dplyr::select(accession, species),
      by = "accession"
    ) |>
    dplyr::full_join(bdqc_grouped, by = "species") |>
    dplyr::filter(!is.na(groupe))),

  tar_target(gene_prevalence_plot, QCGenomLandscape::plot_gene_prevalence(genes_grouped)),
  tar_target(gene_prevalence_saved, {
    ggplot2::ggsave(file.path(results_dir, "genes_prevalence.svg"), gene_prevalence_plot, width = 18, height = 18)
    file.path(results_dir, "genes_prevalence.svg")
  }, format = "file"),

  tar_target(risk_genes, dplyr::bind_rows(ca_risk_en, qc_risk_en) |>
    dplyr::left_join(genes_grouped |> dplyr::distinct(species, gene_group), by = "species") |>
    dplyr::mutate(has_data = !is.na(gene_group))),
  tar_target(risk_status_plot, QCGenomLandscape::plot_risk_status_coverage(risk_genes, jurisdiction = "QC")),
  tar_target(risk_status_saved, {
    ggplot2::ggsave(file.path(results_dir, "risk_status_coverage.svg"), risk_status_plot, width = 16, height = 10)
    file.path(results_dir, "risk_status_coverage.svg")
  }, format = "file"),

  # ---- 6. Summary table ---------------------------------------------------
  tar_target(
    summary_table,
    QCGenomLandscape::build_summary_dataframe(ncbi_voucher$results, gene_annotations, bdqc_taxo, ca_risk, qc_risk)
  ),

  # ---- 7. Sequence QC ------------------------------------------------------
  tar_target(qc_accessions, ncbi_voucher$results |>
    dplyr::filter(!is.na(organism), !is.na(accession)) |>
    dplyr::pull(accession) |>
    unique()),
  tar_target(gb_records, QCGenomLandscape::fetch_gb_records(qc_accessions)),
  tar_target(seq_data, {
    parsed <- purrr::map_dfr(gb_records, QCGenomLandscape::parse_gb_records)
    dplyr::bind_cols(parsed, QCGenomLandscape::score_sequence_quality(parsed$sequence)) |>
      dplyr::mutate(has_stop = purrr::map_lgl(sequence, QCGenomLandscape::has_stop_codon_coi))
  }),
  tar_target(seq_qc_saved, {
    saveRDS(seq_data, file.path(results_dir, "sequence_qc.rds"))
    file.path(results_dir, "sequence_qc.rds")
  }, format = "file"),

  # ---- 8. Arbutus/Swift export ---------------------------------------------
  tar_target(swift_token, get_swift_token()),
  tar_target(bold_swift, maybe_upload_to_swift(bold_saved, swift_token)),
  tar_target(ncbi_results_swift, maybe_upload_to_swift(ncbi_results_saved, swift_token)),
  tar_target(deficient_queries_swift, maybe_upload_to_swift(deficient_queries_saved, swift_token)),
  tar_target(high_id_queries_swift, maybe_upload_to_swift(high_id_queries_saved, swift_token)),
  tar_target(
    ncbi_non_voucher_results_swift,
    maybe_upload_to_swift(ncbi_non_voucher_results_saved, swift_token)
  ),
  tar_target(
    deficient_queries_non_voucher_swift,
    maybe_upload_to_swift(deficient_queries_non_voucher_saved, swift_token)
  ),
  tar_target(
    high_id_queries_non_voucher_swift,
    maybe_upload_to_swift(high_id_queries_non_voucher_saved, swift_token)
  ),
  tar_target(ncbi_genome_results_swift, maybe_upload_to_swift(ncbi_genome_results_saved, swift_token)),
  tar_target(genes_swift, maybe_upload_to_swift(genes_saved, swift_token)),
  tar_target(gene_prevalence_swift, maybe_upload_to_swift(gene_prevalence_saved, swift_token)),
  tar_target(risk_status_swift, maybe_upload_to_swift(risk_status_saved, swift_token)),
  tar_target(seq_qc_swift, maybe_upload_to_swift(seq_qc_saved, swift_token))
)
