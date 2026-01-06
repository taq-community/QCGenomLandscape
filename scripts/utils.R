# Function to parse lat/lon strings and convert to signed decimal degrees
parse_latlon <- function(latlon_string) {
  if (is.na(latlon_string) || latlon_string == "") {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Extract all numbers (including decimals)
  numbers <- stringr::str_extract_all(latlon_string, "\\d+\\.?\\d*")[[1]]

  # Extract cardinal directions
  has_N <- stringr::str_detect(latlon_string, "N")
  has_S <- stringr::str_detect(latlon_string, "S")
  has_E <- stringr::str_detect(latlon_string, "E")
  has_W <- stringr::str_detect(latlon_string, "W")

  if (length(numbers) < 2) {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Parse latitude (first number)
  lat <- as.numeric(numbers[1])
  if (has_S) lat <- -lat  # South is negative

  # Parse longitude (second number)
  lon <- as.numeric(numbers[2])
  if (has_W) lon <- -lon  # West is negative

  return(list(lat = lat, lon = lon))
}

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
