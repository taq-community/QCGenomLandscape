#' The 16 Atlas datasets whose records are already molecular/eDNA-derived
#'
#' Identified by keyword search across all ~880 Atlas dataset names (sequence,
#' DNA, barcode, metagenome, eDNA, INSDC, GenBank, BOLD, ITS, etc.), then
#' manually verified against false positives -- two matches were dropped
#' because the keyword hit was in the *publisher's* institutional name (e.g.
#' "CIBIO (Research Center in Biodiversity and Genetic Resources)"), not the
#' dataset's actual content. Comparing eDNA against the *whole* Atlas would be
#' partly circular, since a meaningful share of it is already eDNA/molecular
#' data -- some of it (iBOL) the same underlying source as this package's own
#' BOLD extraction.
#'
#' @return Character vector of exact `dataset_name` values
#' @export
atlas_molecular_datasets <- function() {
  c(
    "International Barcode of Life project (iBOL)",
    "Trace element-contaminated soil Metagenome",
    "INSDC Sequences",
    "Centre for Biodiversity Genomics - Canadian Specimens",
    "DFO Quebec Gulf of St. Lawrence eDNA metabarcoding (16S) 2020-2022",
    "DFO Quebec Gulf of St. Lawrence eDNA metabarcoding (12S) 2020-2022",
    "NHMO DNA Bank Vascular plants collection",
    "Microbial Community Database (MiCoDa). A curated global 16S rRNA gene amplicon dataset from all environments",
    "INSDC Environment Sample Sequences",
    "INSDC Host Organism Sequences",
    "DFO Quebec Gulf of St. Lawrence eDNA metabarcoding (COI) 2020-2022",
    "Fungal Internal Transcribed Spacer RNA (ITS) RefSeq Targeted Loci Project",
    "University of Tartu Natural History Museum and Botanical Garden  DNA and Environmental Sample Collections",
    "Centre for Biodiversity Genomics (BIOUG) - Marine Invertebrates",
    "NHMO DNA Bank Arthropod collection",
    "Barcoded Reticulariaceae of the World"
  )
}

#' Fetch a Biodiversité Québec Atlas public-data export, caching it locally
#'
#' Downloads via the same public, tokenless S3 bucket documented at
#' `source("http://atlas.biodiversite-quebec.ca/bq-atlas-parquet.R")` --
#' re-implemented directly (rather than `source()`-ing that URL) so the
#' fetch is inspectable and testable. Files are ~2 GB; already-downloaded
#' exports are reused, not re-fetched.
#'
#' @param export_date Character `"YYYY-MM-DD"`, default the latest export
#'   listed in the bucket's `atlas_export_dates.csv`
#' @param dest_dir Directory to cache the downloaded parquet file in,
#'   default `"data/atlas"`
#' @param base_url Base URL of the Atlas parquet exports
#' @param list_dates_fn Function with signature `(url)` returning the export
#'   dates as a character vector, default reads `atlas_export_dates.csv`;
#'   injectable for testing without network access
#' @param download_fn Function with signature `(url, destfile)`, default
#'   [utils::download.file()]; injectable for testing without network access
#' @return Character scalar, local path to the parquet file
#' @export
fetch_atlas_parquet <- function(export_date = NULL,
                                 dest_dir = "data/atlas",
                                 base_url = "https://object-arbutus.cloud.computecanada.ca/bq-io/atlas/parquet/",
                                 list_dates_fn = function(url) {
                                   dates <- utils::read.csv(url, header = FALSE, col.names = "date")$date
                                   dates[nzchar(trimws(dates))]
                                 },
                                 download_fn = function(url, destfile) {
                                   utils::download.file(url, destfile, mode = "wb", quiet = FALSE)
                                 }) {
  if (is.null(export_date)) {
    dates <- list_dates_fn(paste0(base_url, "atlas_export_dates.csv"))
    export_date <- max(dates)
    logger::log_info("fetch_atlas_parquet: latest export date is {export_date}")
  }

  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  file_name <- sprintf("atlas_public_%s.parquet", export_date)
  dest <- file.path(dest_dir, file_name)

  if (!file.exists(dest)) {
    logger::log_info("fetch_atlas_parquet: downloading {file_name} (~2 GB)...")
    download_fn(paste0(base_url, file_name), dest)
    logger::log_success("fetch_atlas_parquet: saved to {dest}")
  } else {
    logger::log_info("fetch_atlas_parquet: using cached {dest}")
  }

  dest
}

#' Build a full hexagonal grid over a boundary, with stable cell IDs
#'
#' Unlike [make_hex_grid()] (which counts one point layer and drops empty
#' cells), this returns *every* cell with a stable `cell_id` regardless of
#' occupancy, so multiple point layers (e.g. eDNA vs traditional
#' occurrences) can be joined against the same reusable grid and compared
#' cell-for-cell.
#'
#' @param boundary_sf `sf` polygon object defining the grid extent
#' @param cellsize Numeric, grid cell size in meters (in `crs_proj`), default 10000
#' @param crs_proj Numeric EPSG code for the projection used to build the
#'   grid, default `3347` (Canada Lambert Conformal Conic)
#' @return `sf` polygon grid (WGS84) with a `cell_id` column
#' @export
build_hex_grid_cells <- function(boundary_sf, cellsize = 10000, crs_proj = 3347) {
  boundary_proj <- sf::st_transform(boundary_sf, crs_proj)

  hex_grid <- sf::st_make_grid(boundary_proj, cellsize = cellsize, square = FALSE) |>
    sf::st_as_sf() |>
    sf::st_intersection(sf::st_union(boundary_proj))
  hex_grid$cell_id <- seq_len(nrow(hex_grid))

  sf::st_transform(hex_grid, 4326)
}

#' Compare eDNA vs traditional occurrence coverage against the Atlas, on a hex grid
#'
#' The heavy lifting: classifies each Atlas record as `"eDNA"` (one of
#' [atlas_molecular_datasets()]) or `"traditional"`, spatially bins every
#' record into `hex_grid` via a DuckDB spatial join (fast even at Atlas
#' scale -- tens of millions of rows join against a ~18k-cell grid in well
#' under a minute), and folds in `edna_occurrences` (e.g. from
#' [extract_edna_occurrences()]) as additional eDNA evidence, binned into
#' the same grid via [sf::st_join()].
#'
#' No deduplication is attempted between `edna_occurrences` and Atlas's own
#' molecular datasets (some of which, e.g. iBOL, share underlying records
#' with a BOLD-sourced `edna_occurrences`) -- this has little effect on the
#' presence-based metrics [summarize_edna_contribution()] computes from the
#' result, which count distinct (species, cell) / (species, year) pairs, not
#' raw record volume.
#'
#' @param atlas_parquet_path Character, local path to an Atlas parquet
#'   export, e.g. from [fetch_atlas_parquet()]
#' @param hex_grid `sf` polygon grid with a `cell_id` column, e.g. from
#'   [build_hex_grid_cells()]
#' @param edna_occurrences Tibble with columns `species`, `lon`, `lat`,
#'   `date`, e.g. from [extract_edna_occurrences()]
#' @param bdqc_species Character vector of species to restrict the
#'   comparison to, or `NULL` for no restriction (default `NULL`)
#' @param molecular_datasets Character vector of Atlas `dataset_name` values
#'   to treat as eDNA, default [atlas_molecular_datasets()]
#' @param con An existing DBI connection to a DuckDB database, or `NULL`
#'   (default) to open and close a temporary one; injectable for testing
#' @return A list with `edna_agg` and `traditional_agg` tibbles, each with
#'   columns `species`, `cell_id`, `year_obs`, `n`
#' @details Requires the DuckDB `httpfs` and `spatial` extensions. Install them
#'   once with `duckdb::duckdb_install_extension(c("httpfs", "spatial"))` before
#'   first use.
#' @export
compare_edna_atlas_coverage <- function(atlas_parquet_path,
                                         hex_grid,
                                         edna_occurrences,
                                         bdqc_species = NULL,
                                         molecular_datasets = atlas_molecular_datasets(),
                                         con = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
    stop("compare_edna_atlas_coverage() needs the DBI and duckdb packages -- install.packages(c(\"DBI\", \"duckdb\"))", call. = FALSE)
  }

  own_con <- is.null(con)
  if (own_con) {
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }
  DBI::dbExecute(con, "LOAD httpfs;")
  DBI::dbExecute(con, "LOAD spatial;")

  hex_gpkg <- tempfile(fileext = ".gpkg")
  on.exit(unlink(hex_gpkg), add = TRUE)
  sf::st_write(hex_grid["cell_id"], hex_gpkg, delete_dsn = TRUE, quiet = TRUE)
  DBI::dbExecute(con, sprintf("CREATE OR REPLACE TABLE hexgrid AS SELECT * FROM ST_Read('%s');", hex_gpkg))
  DBI::dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW atlas AS SELECT * FROM read_parquet('%s');",
    atlas_parquet_path
  ))

  in_clause <- paste(sprintf("'%s'", gsub("'", "''", molecular_datasets)), collapse = ", ")
  species_filter <- ""
  if (!is.null(bdqc_species)) {
    species_in <- paste(sprintf("'%s'", gsub("'", "''", bdqc_species)), collapse = ", ")
    species_filter <- sprintf("AND a.valid_scientific_name IN (%s)", species_in)
  }

  logger::log_info("compare_edna_atlas_coverage: spatial join + aggregation against Atlas parquet...")
  atlas_agg <- DBI::dbGetQuery(con, sprintf("
    SELECT
      a.valid_scientific_name AS species,
      CASE WHEN a.dataset_name IN (%s) THEN 'eDNA' ELSE 'traditional' END AS group_type,
      h.cell_id,
      a.year_obs,
      count(*) AS n
    FROM atlas a
    JOIN hexgrid h ON ST_Intersects(ST_Point(a.longitude, a.latitude), h.geom)
    WHERE a.longitude IS NOT NULL AND a.latitude IS NOT NULL AND a.valid_scientific_name IS NOT NULL
      %s
    GROUP BY a.valid_scientific_name, group_type, h.cell_id, a.year_obs
  ", in_clause, species_filter)) |> tibble::as_tibble()
  logger::log_success("compare_edna_atlas_coverage: {nrow(atlas_agg)} aggregated rows from Atlas")

  pts <- sf::st_as_sf(edna_occurrences, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  joined <- suppressWarnings(sf::st_join(pts, hex_grid["cell_id"], join = sf::st_intersects))
  own_agg <- joined |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(.data$cell_id)) |>
    dplyr::mutate(year_obs = as.integer(format(.data$date, "%Y"))) |>
    dplyr::count(.data$species, .data$cell_id, .data$year_obs, name = "n")

  if (!is.null(bdqc_species)) {
    atlas_agg <- atlas_agg |> dplyr::filter(.data$species %in% bdqc_species)
    own_agg <- own_agg |> dplyr::filter(.data$species %in% bdqc_species)
  }

  edna_agg <- dplyr::bind_rows(
    atlas_agg |> dplyr::filter(.data$group_type == "eDNA") |> dplyr::select("species", "cell_id", "year_obs", "n"),
    own_agg
  ) |>
    dplyr::group_by(.data$species, .data$cell_id, .data$year_obs) |>
    dplyr::summarise(n = sum(.data$n), .groups = "drop")

  traditional_agg <- atlas_agg |>
    dplyr::filter(.data$group_type == "traditional") |>
    dplyr::select("species", "cell_id", "year_obs", "n")

  list(edna_agg = edna_agg, traditional_agg = traditional_agg)
}

#' Summarize an eDNA-vs-traditional coverage comparison into contribution metrics
#'
#' Pure post-processing over the output of [compare_edna_atlas_coverage()] --
#' no I/O, so it's cheap to recompute (e.g. per taxonomic subgroup) without
#' re-running the spatial join.
#'
#' @param edna_agg,traditional_agg Tibbles with columns `species`, `cell_id`,
#'   `year_obs`, `n`, e.g. from [compare_edna_atlas_coverage()]
#' @return A list: `edna_species`, `traditional_species`, `edna_only_species`,
#'   `traditional_only_species`, `shared_species` (character vectors);
#'   `novel_cells`, `novel_years` (tibbles of `species` + the novel
#'   cell/year); and summary counts/percentages (`n_species_edna`,
#'   `n_species_edna_only`, `pct_species_edna_only`, `n_species_traditional`,
#'   `n_species_traditional_only`, `pct_species_traditional_only`,
#'   `n_species_shared`, `jaccard_species` (0 = disjoint species sets, 1 =
#'   identical), `n_cells_edna`, `n_cells_novel`, `pct_cells_novel`,
#'   `n_years_edna`, `n_years_novel`, `pct_years_novel`)
#' @export
summarize_edna_contribution <- function(edna_agg, traditional_agg) {
  edna_species <- unique(edna_agg$species)
  traditional_species <- unique(traditional_agg$species)
  edna_only_species <- setdiff(edna_species, traditional_species)
  traditional_only_species <- setdiff(traditional_species, edna_species)
  shared_species <- intersect(edna_species, traditional_species)
  all_species <- union(edna_species, traditional_species)

  edna_cells <- dplyr::distinct(edna_agg, .data$species, .data$cell_id)
  traditional_cells <- dplyr::distinct(traditional_agg, .data$species, .data$cell_id)
  novel_cells <- dplyr::anti_join(edna_cells, traditional_cells, by = c("species", "cell_id"))

  edna_years <- edna_agg |> dplyr::filter(!is.na(.data$year_obs)) |> dplyr::distinct(.data$species, .data$year_obs)
  traditional_years <- traditional_agg |>
    dplyr::filter(!is.na(.data$year_obs)) |>
    dplyr::distinct(.data$species, .data$year_obs)
  novel_years <- dplyr::anti_join(edna_years, traditional_years, by = c("species", "year_obs"))

  list(
    edna_species = edna_species,
    traditional_species = traditional_species,
    edna_only_species = edna_only_species,
    traditional_only_species = traditional_only_species,
    shared_species = shared_species,
    novel_cells = novel_cells,
    novel_years = novel_years,
    n_species_edna = length(edna_species),
    n_species_edna_only = length(edna_only_species),
    pct_species_edna_only = 100 * length(edna_only_species) / length(edna_species),
    n_species_traditional = length(traditional_species),
    n_species_traditional_only = length(traditional_only_species),
    pct_species_traditional_only = 100 * length(traditional_only_species) / length(traditional_species),
    n_species_shared = length(shared_species),
    # indice de Jaccard sur les ensembles d'especes -- 0 = ensembles disjoints
    # (diversite totalement differente), 1 = ensembles identiques
    jaccard_species = length(shared_species) / length(all_species),
    n_cells_edna = nrow(edna_cells),
    n_cells_novel = nrow(novel_cells),
    pct_cells_novel = 100 * nrow(novel_cells) / nrow(edna_cells),
    n_years_edna = nrow(edna_years),
    n_years_novel = nrow(novel_years),
    pct_years_novel = 100 * nrow(novel_years) / nrow(edna_years)
  )
}
