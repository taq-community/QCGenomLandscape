#' Load a CanVec administrative-boundary layer, optionally filtered
#'
#' @param gdb_path Character, path to the `.gdb` directory (e.g.
#'   `"data/canvec_1M_CA_Admin.gdb"`)
#' @param layer Character, layer name, default `"geo_political_region_2"`
#' @param jurisdiction Integer or `NULL`; if set, keeps rows where
#'   `jurisdiction == jurisdiction` (e.g. `102` for Quebec)
#' @param country Integer or `NULL`; if set, keeps rows where
#'   `country == country` (e.g. `140` for Canada)
#' @param crs Target CRS, default `4326`
#' @return An `sf` object
#' @export
load_canvec_boundary <- function(gdb_path, layer = "geo_political_region_2",
                                  jurisdiction = NULL, country = NULL, crs = 4326) {
  boundary <- sf::read_sf(gdb_path, layer = layer)

  if (!is.null(jurisdiction)) {
    boundary <- dplyr::filter(boundary, jurisdiction == .env$jurisdiction)
  }
  if (!is.null(country)) {
    boundary <- dplyr::filter(boundary, country == .env$country)
  }

  sf::st_transform(boundary, crs)
}

#' Flag points as within a polygon boundary
#'
#' @param points_sf `sf` points object
#' @param boundary_sf `sf` polygon object (transformed to `points_sf`'s CRS
#'   if it differs)
#' @param flag_name Character, name of the new logical column, default `"in_boundary"`
#' @return `points_sf` with an added logical column named `flag_name`
#' @export
flag_within_boundary <- function(points_sf, boundary_sf, flag_name = "in_boundary") {
  if (sf::st_crs(boundary_sf) != sf::st_crs(points_sf)) {
    boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(points_sf))
  }
  points_sf[[flag_name]] <- lengths(sf::st_within(points_sf, boundary_sf)) > 0
  points_sf
}

#' Normalize NCBI + BOLD records into a common eDNA occurrence table
#'
#' Pulls the fields needed to compare eDNA records against traditional
#' occurrence data (e.g. the Biodiversité Québec Atlas) spatio-temporally:
#' one row per geo-referenced eDNA record, with species, coordinates, and
#' collection date normalized to a common shape regardless of source.
#' Records with no parseable coordinates, or no species name, are dropped --
#' both sources have some (BOLD: unidentified specimens; NCBI: missing
#' `lat_lon`).
#'
#' @param ncbi_results Tibble as returned by [fetch_ncbi_sequences()]'s
#'   `results` (must have `organism`, `lat_lon`, `collection_date`,
#'   `accession` columns)
#' @param bold_raw Character scalar, raw BOLD TSV as returned by
#'   [fetch_bold_sequences()]
#' @return Tibble with columns `source` (`"NCBI"`/`"BOLD"`), `species`,
#'   `lon`, `lat`, `date` (`Date`, `NA` where uncollectable), `record_id`
#' @export
extract_edna_occurrences <- function(ncbi_results, bold_raw) {
  ncbi_occ <- ncbi_results |>
    dplyr::filter(!is.na(organism), !is.na(lat_lon)) |>
    dplyr::mutate(parsed = purrr::map(lat_lon, parse_latlon)) |>
    dplyr::transmute(
      source = "NCBI",
      species = organism,
      lon = purrr::map_dbl(parsed, "lon"),
      lat = purrr::map_dbl(parsed, "lat"),
      date = parse_gb_collection_date(collection_date),
      record_id = accession
    ) |>
    dplyr::filter(!is.na(lon), !is.na(lat))

  bold_df <- readr::read_tsv(
    I(bold_raw),
    col_types = readr::cols_only(
      species = readr::col_character(),
      coord = readr::col_character(),
      collection_date_start = readr::col_character(),
      record_id = readr::col_character()
    ),
    progress = FALSE
  )
  bold_occ <- bold_df |>
    dplyr::filter(!is.na(species), species != "", !is.na(coord)) |>
    dplyr::mutate(parsed = purrr::map(coord, parse_bold_coord)) |>
    dplyr::transmute(
      source = "BOLD",
      species,
      lon = purrr::map_dbl(parsed, "lon"),
      lat = purrr::map_dbl(parsed, "lat"),
      date = parse_gb_collection_date(collection_date_start),
      record_id = record_id
    ) |>
    dplyr::filter(!is.na(lon), !is.na(lat))

  dplyr::bind_rows(ncbi_occ, bold_occ)
}

#' Build a hexagonal summary grid of point counts within a boundary
#'
#' @param points_sf `sf` points object (any CRS)
#' @param boundary_sf `sf` polygon object defining the grid extent
#' @param cellsize Numeric, grid cell size in meters (in `crs_proj`), default 10000
#' @param crs_proj Numeric EPSG code for the projection used to build the
#'   grid, default `3347` (Canada Lambert Conformal Conic)
#' @param log_scale Logical; if `TRUE`, adds a `log_count = log(n_records + 1)` column
#' @return `sf` polygon grid (WGS84), filtered to cells with `n_records > 0`
#' @export
make_hex_grid <- function(points_sf, boundary_sf, cellsize = 10000,
                           crs_proj = 3347, log_scale = TRUE) {
  points_proj <- sf::st_transform(points_sf, crs_proj)
  boundary_proj <- sf::st_transform(boundary_sf, crs_proj)

  hex_grid <- sf::st_make_grid(boundary_proj, cellsize = cellsize, square = FALSE) |>
    sf::st_as_sf() |>
    sf::st_intersection(boundary_proj)

  hex_grid$n_records <- lengths(sf::st_intersects(hex_grid, points_proj))
  hex_grid <- hex_grid[hex_grid$n_records > 0, ]

  if (log_scale) {
    hex_grid$log_count <- log(hex_grid$n_records + 1)
  }

  sf::st_transform(hex_grid, 4326)
}
