# Utility functions for QCGenomLandscape
# -----------------------------------------------------------------------------

#' Parse lat/lon strings from NCBI format to decimal degrees
#'
#' @param latlon_string A string like "45.5 N 73.6 W" or "45.5N 73.6W"
#' @return A list with lat and lon as numeric values
parse_latlon <- function(latlon_string) {
  if (is.na(latlon_string) || latlon_string == "") {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Extract all numbers (including decimals)
  numbers <- stringr::str_extract_all(latlon_string, "\\d+\\.?\\d*")[[1]]

  # Extract cardinal directions
  has_S <- stringr::str_detect(latlon_string, "S")
  has_W <- stringr::str_detect(latlon_string, "W")

  if (length(numbers) < 2) {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Parse latitude (first number)
  lat <- as.numeric(numbers[1])
  if (has_S) lat <- -lat # South is negative

  # Parse longitude (second number)
  lon <- as.numeric(numbers[2])
  if (has_W) lon <- -lon # West is negative

  return(list(lat = lat, lon = lon))
}
