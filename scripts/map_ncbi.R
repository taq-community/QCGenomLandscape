source("scripts/utils.R")

# Parse lat/lon and add as separate columns
ncbi_results <- readRDS("results/ncbi_results.rds") |>
  dplyr::mutate(
    parsed_coords = purrr::map(lat_lon, parse_latlon),
    latitude = purrr::map_dbl(parsed_coords, "lat"),
    longitude = purrr::map_dbl(parsed_coords, "lon")
  ) |>
  dplyr::select(-parsed_coords)

ncbi_sf <- sf::st_as_sf(
  dplyr::filter(ncbi_results, !is.na(latitude) & !is.na(longitude)),
  coords = c("longitude", "latitude"),
  crs = 4326
)

# From https://open.canada.ca/data/en/dataset/306e5004-534b-4110-9feb-58e3a5c3fd97
qc <- sf::read_sf("data/canvec_1M_CA_Admin.gdb", layer = "geo_political_region_2") |>
  dplyr::filter(jurisdiction == 102) |>
  sf::st_transform(4326)

can <- sf::read_sf("data/canvec_1M_CA_Admin.gdb", layer = "geo_political_region_2") |>
  dplyr::filter(country == 140) |>
  sf::st_transform(4326)

ncbi_sf <- ncbi_sf |> dplyr::mutate(
  in_ca = lengths(sf::st_within(geometry, can)) > 0,
  in_qc = lengths(sf::st_within(geometry, qc)) > 0
)

mapview::mapview(ncbi_sf)
mapview::mapview(ncbi_sf |> dplyr::filter(in_ca | in_qc))

# Create hexagon grid at 10km resolution for Canada
ncbi_ca <- ncbi_sf |> dplyr::filter(in_qc)

# Transform to a projected CRS for accurate distance calculations (Canada Lambert Conformal Conic)
ncbi_ca_proj <- sf::st_transform(ncbi_ca, 3347)
can_proj <- sf::st_transform(can, 3347)

# Create hexagonal grid over Canada (10km = 10000m cell size)
hex_grid <- sf::st_make_grid(
  can_proj,
  cellsize = 50000,
  square = FALSE
) |>
  sf::st_as_sf() |>
  sf::st_intersection(can_proj)

# Count points in each hexagon
hex_grid$n_records <- lengths(sf::st_intersects(hex_grid, ncbi_ca_proj))

# Filter to hexagons with at least one record
hex_with_data <- hex_grid |>
  dplyr::filter(n_records > 0) 

# Transform back to WGS84 for mapping
hex_with_data <- sf::st_transform(hex_with_data, 4326)

mapview::mapview(hex_with_data, zcol = "n_records")
