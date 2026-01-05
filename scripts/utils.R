
# Function to parse lat/lon strings and convert to signed decimal degrees
parse_latlon <- function(latlon_string) {
  if (is.na(latlon_string) || latlon_string == "") {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Extract all numbers (including decimals)
  numbers <- str_extract_all(latlon_string, "\\d+\\.?\\d*")[[1]]

  # Extract cardinal directions
  has_N <- str_detect(latlon_string, "N")
  has_S <- str_detect(latlon_string, "S")
  has_E <- str_detect(latlon_string, "E")
  has_W <- str_detect(latlon_string, "W")

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
ncbi_results <- ncbi_results |>
  mutate(
    parsed_coords = map(lat_lon, parse_latlon),
    latitude = map_dbl(parsed_coords, "lat"),
    longitude = map_dbl(parsed_coords, "lon")
  ) |>
  select(-parsed_coords)

ncbi_sf <- sf::st_as_sf(
  dplyr::filter(ncbi_results, !is.na(latitude) & !is.na(longitude)),
  coords = c("longitude", "latitude"), 
  crs = 4326
)

mapview::mapview(ncbi_sf)

### GENOME AVAILABLE ###
# Vérifier s'il existe des génomes complets
# genome_search <- rentrez::entrez_search(
#   db = "genome",
#   term = species_name,
#   retmax = 0
# )
