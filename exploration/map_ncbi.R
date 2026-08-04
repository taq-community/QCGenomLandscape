# Quick interactive look at NCBI Quebec/Canada records + a hex-density map.
# Not part of the reproducible _targets.R pipeline.
library(QCGenomLandscape)

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
qc <- load_canvec_boundary("data/canvec_1M_CA_Admin.gdb", jurisdiction = 102)
can <- load_canvec_boundary("data/canvec_1M_CA_Admin.gdb", country = 140)

ncbi_sf <- ncbi_sf |>
  flag_within_boundary(can, "in_ca") |>
  flag_within_boundary(qc, "in_qc")

mapview::mapview(ncbi_sf)
mapview::mapview(ncbi_sf |> dplyr::filter(in_ca | in_qc))

# 10km hex-density grid over Quebec
ncbi_qc <- ncbi_sf |> dplyr::filter(in_qc)
hex_with_data <- make_hex_grid(ncbi_qc, qc, cellsize = 10000)

mapview::mapview(hex_with_data, zcol = "n_records")
