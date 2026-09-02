library(targets)
library(tarchetypes)
library(QCGenomLandscape)

tar_option_set(
  packages = c("QCGenomLandscape", "sf", "Biostrings"),
  # packages = just attaches namespaces; imports = is the separate option
  # that actually tracks QCGenomLandscape's function bodies for invalidation
  # -- without it, targets never notices when R/ code changes and silently
  # keeps stale cached results (how ncbi_queries/ncbi_sequences stayed on
  # the pre-merge voucher-only query for weeks after the code changed).
  imports = "QCGenomLandscape",
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

# ---- Cloud-backed store for expensive-to-recompute targets ----------------
# ncbi_sequences/gb_records/seq_data are the multi-hour steps -- storing
# their data in Arbutus S3 (not just local _targets/objects/) lets a
# collaborator fetch already-computed results via tar_make() instead of
# re-running the NCBI fetch/parse from scratch.
#
# Needs AWS_ACCESS_KEY_ID/AWS_SECRET_ACCESS_KEY -- Arbutus S3 "EC2
# credentials", generated with `openstack ec2 credentials create` (or via
# Horizon: Project > API Access > EC2 Credentials), NOT the
# SWIFT_APP_CRED_ID/SECRET used for the native Swift upload step below --
# see https://docs.alliancecan.ca/wiki/Accessing_object_storage_with_s3cmd.
# QCGL_S3_BUCKET must already exist (`s3cmd mb s3://BUCKET_NAME/`).
#
# Falls back to the default local-only store when these aren't set, so
# collaborators without S3 access can still run the pipeline.
use_s3_store <- all(nzchar(Sys.getenv(c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY"))))
s3_repository <- if (use_s3_store) "aws" else "local"
s3_resources <- if (use_s3_store) {
  tar_resources(aws = tar_resources_aws(
    bucket = Sys.getenv("QCGL_S3_BUCKET", "qcgenomlandscape"),
    prefix = "targets",
    region = "us-east-1",
    endpoint = "https://object-arbutus.alliancecan.ca",
    # Arbutus's S3 gateway uses path-style bucket addressing
    # (object-arbutus.alliancecan.ca/BUCKET/KEY), not virtual-hosted-style
    # (BUCKET.object-arbutus.alliancecan.ca/KEY)
    s3_force_path_style = TRUE
  ))
} else {
  logger::log_info("AWS_ACCESS_KEY_ID/AWS_SECRET_ACCESS_KEY not set -- expensive targets use local storage only")
  tar_resources()
}

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
  # A failed upload (expired token, transient network issue, Arbutus outage,
  # ...) shouldn't crash a multi-hour tar_make() run over what's meant to be
  # an optional export step -- same non-fatal treatment as a missing token.
  # upload_to_swift() already logs the error before this catches it.
  tryCatch(
    upload_to_swift(path, token = token),
    error = function(e) {
      warning("Skipping Swift upload of ", path, " after failure: ", e$message, call. = FALSE)
      invisible(NULL)
    }
  )
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

  # ---- 2. NCBI sequences + geo-filter --------------------------------------
  # batch_size batches species sharing a marker into one OR'd Entrez query --
  # cuts total NCBI requests by roughly that factor. Species are recovered
  # afterwards from the `organism` field NCBI returns, not by parsing `query`.
  # No more server-side voucher[Title] restriction (that used to mean two
  # full passes, filtered and unfiltered) -- fetch_ncbi_sequences() infers
  # `is_voucher` per record from its title instead, so one pass covers both.
  tar_target(
    ncbi_queries,
    QCGenomLandscape::build_ncbi_queries(bdqc_species, query_primers, batch_size = 25)
  ),
  tar_target(
    ncbi_sequences,
    QCGenomLandscape::fetch_ncbi_sequences(ncbi_queries),
    repository = s3_repository, resources = s3_resources
  ),
  tar_target(ncbi_results_saved, {
    saveRDS(ncbi_sequences$results, file.path(results_dir, "ncbi_results.rds"))
    file.path(results_dir, "ncbi_results.rds")
  }, format = "file"),
  tar_target(deficient_queries_saved, {
    saveRDS(ncbi_sequences$deficient_queries, file.path(results_dir, "deficient_queries.rds"))
    file.path(results_dir, "deficient_queries.rds")
  }, format = "file"),
  tar_target(high_id_queries_saved, {
    saveRDS(ncbi_sequences$high_id_queries, file.path(results_dir, "high_id_queries.rds"))
    file.path(results_dir, "high_id_queries.rds")
  }, format = "file"),
  # Geo-filtered (QC/CA) view -- demonstrates load_canvec_boundary()/
  # flag_within_boundary(), not saved to a named results/ file since the
  # original script only used this interactively
  tar_target(ncbi_geo, {
    pts <- ncbi_sequences$results |>
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

  # ---- 3. NCBI full genomes -----------------------------------------------
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

  # ---- 4. Gene annotations (subsample) + taxon_representation figures ----
  tar_target(gene_accessions, {
    set.seed(42)
    ncbi_sequences$results |>
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
      ncbi_sequences$results |>
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

  # ---- 5. Summary table ---------------------------------------------------
  tar_target(
    summary_table,
    QCGenomLandscape::build_summary_dataframe(ncbi_sequences$results, gene_annotations, bdqc_taxo, ca_risk, qc_risk)
  ),

  # ---- 6. Sequence QC ------------------------------------------------------
  tar_target(qc_accessions, ncbi_sequences$results |>
    dplyr::filter(!is.na(organism), !is.na(accession)) |>
    dplyr::pull(accession) |>
    unique()),
  tar_target(
    gb_records,
    QCGenomLandscape::fetch_gb_records(qc_accessions),
    repository = s3_repository, resources = s3_resources
  ),
  tar_target(
    seq_data,
    QCGenomLandscape::build_sequence_qc_table(gb_records),
    repository = s3_repository, resources = s3_resources
  ),
  tar_target(seq_qc_saved, {
    saveRDS(seq_data, file.path(results_dir, "sequence_qc.rds"))
    file.path(results_dir, "sequence_qc.rds")
  }, format = "file"),

  # ---- 7. Arbutus/Swift export ---------------------------------------------
  # cue = "always": get_swift_token() has no target-level dependencies, so
  # targets has no way to notice a cached token has simply expired -- force
  # a fresh Keystone auth every run instead of reusing a stale one
  tar_target(swift_token, get_swift_token(), cue = tar_cue(mode = "always")),
  tar_target(bold_swift, maybe_upload_to_swift(bold_saved, swift_token)),
  tar_target(ncbi_results_swift, maybe_upload_to_swift(ncbi_results_saved, swift_token)),
  tar_target(deficient_queries_swift, maybe_upload_to_swift(deficient_queries_saved, swift_token)),
  tar_target(high_id_queries_swift, maybe_upload_to_swift(high_id_queries_saved, swift_token)),
  tar_target(ncbi_genome_results_swift, maybe_upload_to_swift(ncbi_genome_results_saved, swift_token)),
  tar_target(genes_swift, maybe_upload_to_swift(genes_saved, swift_token)),
  tar_target(gene_prevalence_swift, maybe_upload_to_swift(gene_prevalence_saved, swift_token)),
  tar_target(risk_status_swift, maybe_upload_to_swift(risk_status_saved, swift_token)),
  tar_target(seq_qc_swift, maybe_upload_to_swift(seq_qc_saved, swift_token)),

  # ---- 8. Marker data description ------------------------------------------
  # Extract the query_marker string (everything after " AND " in the Entrez
  # query) to identify which marker group each record belongs to, then join
  # back to query_primers for human-readable labels.
  #
  # Some groups share the same query_marker (Birds/Fish/Mammals/Reptiles all
  # use the same COI+12S+16S+cytb query; Conifers/Other plants share the
  # rbcL+matK+ITS+trnL query). Pre-aggregate to one row per query_marker so
  # the left_join doesn't fan out into many-to-many duplicates.
  tar_target(marker_stats, {
    primers_agg <- query_primers |>
      dplyr::group_by(query_marker) |>
      dplyr::summarise(
        groups  = paste(sort(unique(group)), collapse = " / "),
        markers = dplyr::first(markers),
        .groups = "drop"
      )

    ncbi_sequences$results |>
      dplyr::mutate(
        query_marker = stringr::str_extract(query, "(?<= AND ).+$")
      ) |>
      dplyr::left_join(primers_agg, by = "query_marker") |>
      dplyr::group_by(groups, markers, query_marker) |>
      dplyr::summarise(
        n_sequences = dplyr::n(),
        n_species    = dplyr::n_distinct(organism, na.rm = TRUE),
        median_slen  = stats::median(slen, na.rm = TRUE),
        sd_slen      = stats::sd(slen, na.rm = TRUE),
        q25_slen     = stats::quantile(slen, 0.25, na.rm = TRUE),
        q75_slen     = stats::quantile(slen, 0.75, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(n_sequences))
  }),
  tar_target(marker_stats_saved, {
    saveRDS(marker_stats, file.path(results_dir, "marker_stats.rds"))
    file.path(results_dir, "marker_stats.rds")
  }, format = "file"),

  # ---- 9. Intraspecific sequence-length variation per marker ---------------
  # For each (species, marker) with >1 record, compute length variability.
  # High CV (sd / median * 100) signals mixed-length sequences for the same
  # species+marker — likely lower-quality submissions or mixed amplicons.
  #
  # slen <= 10000 excludes complete genomes/mitogenomes (15kb+) that NCBI
  # returns for gene queries (COI[Gene] also matches whole mitogenomes); those
  # inflate CV to meaningless values and aren't barcode-quality candidates.
  tar_target(intraspecific_variation, {
    primers_dedup <- dplyr::distinct(query_primers, query_marker, markers)

    ncbi_sequences$results |>
      dplyr::filter(!is.na(organism), !is.na(slen), slen <= 10000) |>
      dplyr::mutate(
        query_marker = stringr::str_extract(query, "(?<= AND ).+$")
      ) |>
      dplyr::left_join(primers_dedup, by = "query_marker") |>
      dplyr::group_by(organism, markers) |>
      dplyr::filter(dplyr::n() > 1) |>
      dplyr::summarise(
        n_seq       = dplyr::n(),
        median_slen = stats::median(slen),
        sd_slen     = stats::sd(slen),
        cv_pct      = round(sd_slen / median_slen * 100, 1),
        min_slen    = min(slen),
        max_slen    = max(slen),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(cv_pct))
  }),
  tar_target(intraspecific_variation_saved, {
    saveRDS(intraspecific_variation, file.path(results_dir, "intraspecific_variation.rds"))
    file.path(results_dir, "intraspecific_variation.rds")
  }, format = "file")
)
