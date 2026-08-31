#' Parse an NCBI-style lat/lon string into decimal degrees
#'
#' @param latlon_string A string like "45.5 N 73.6 W" or "45.5N 73.6W"
#' @return A list with `lat` and `lon` as numeric values (`NA_real_` if the
#'   string is missing, empty, or doesn't contain at least two numbers)
#' @examples
#' parse_latlon("45.5 N 73.6 W")
#' parse_latlon("33.9 S 18.4 E")
#' parse_latlon(NA)
#' @export
parse_latlon <- function(latlon_string) {
  if (is.na(latlon_string) || latlon_string == "") {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Extract all numbers (including decimals)
  numbers <- regmatches(latlon_string, gregexpr("\\d+\\.?\\d*", latlon_string))[[1]]

  # Extract cardinal directions
  has_S <- grepl("S", latlon_string, fixed = TRUE)
  has_W <- grepl("W", latlon_string, fixed = TRUE)

  if (length(numbers) < 2) {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  # Parse latitude (first number)
  lat <- as.numeric(numbers[1])
  if (has_S) lat <- -lat # South is negative

  # Parse longitude (second number)
  lon <- as.numeric(numbers[2])
  if (has_W) lon <- -lon # West is negative

  list(lat = lat, lon = lon)
}

#' Parse a BOLD-style `coord` string into decimal degrees
#'
#' @param coord_string A string like `"[45.45, -75.77]"` (`"[lat, lon]"`)
#' @return A list with `lat` and `lon` as numeric values (`NA_real_` if the
#'   string is missing, empty, or doesn't contain exactly two numbers)
#' @examples
#' parse_bold_coord("[45.45, -75.77]")
#' parse_bold_coord(NA)
#' @export
parse_bold_coord <- function(coord_string) {
  if (is.na(coord_string) || coord_string == "") {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  numbers <- regmatches(coord_string, gregexpr("-?\\d+\\.?\\d*", coord_string))[[1]]
  if (length(numbers) != 2) {
    return(list(lat = NA_real_, lon = NA_real_))
  }

  list(lat = as.numeric(numbers[1]), lon = as.numeric(numbers[2]))
}

#' Parse a GenBank/BOLD-style collection-date string into a `Date`
#'
#' Collection dates in GenBank (and, less often, BOLD) records come in
#' several formats of varying precision -- `"DD-Mon-YYYY"` (the common
#' case), `"Mon-YYYY"`/`"YYYY-MM"`/`"YYYY"` (partial dates, normalized to
#' the first of the month/year), plus ISO `"YYYY-MM-DD"`. Occasionally the
#' field holds non-date junk instead (e.g. a locality string or a
#' coordinate pair that leaked into it upstream); that returns `NA`.
#'
#' `as.Date()` is not used directly here: with an explicit `format=` it
#' does *partial*, not full-string, matching (`as.Date("58.298 N 8.539 E",
#' format = "%Y")` silently parses to the year 58 using today's month/day)
#' and formats lacking `%d` (`"%b-%Y"`, `"%Y-%m"`) fail outright rather
#' than defaulting the day. Each recognized pattern is instead converted to
#' an explicit ISO string first, so the final parse is unambiguous.
#'
#' @param x Character vector of collection-date values
#' @return `Date` vector, the same length as `x`
#' @examples
#' parse_gb_collection_date(c("28-Sep-2014", "Sep-2023", "2020", "junk", NA))
#' @export
parse_gb_collection_date <- function(x) {
  iso <- vapply(x, function(val) {
    if (is.na(val)) {
      return(NA_character_)
    }
    val <- trimws(val)

    if (grepl("^\\d{2}-[A-Za-z]{3}-\\d{4}$", val)) {
      d <- as.Date(val, format = "%d-%b-%Y")
      return(if (is.na(d)) NA_character_ else format(d, "%Y-%m-%d"))
    }
    if (grepl("^\\d{4}-\\d{2}-\\d{2}$", val)) {
      return(val)
    }
    if (grepl("^\\d{4}-\\d{2}$", val)) {
      return(paste0(val, "-01"))
    }
    if (grepl("^[A-Za-z]{3}-\\d{4}$", val)) {
      d <- as.Date(paste0("01-", val), format = "%d-%b-%Y")
      return(if (is.na(d)) NA_character_ else format(d, "%Y-%m-%d"))
    }
    if (grepl("^\\d{4}$", val)) {
      return(paste0(val, "-01-01"))
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)

  as.Date(iso, format = "%Y-%m-%d")
}
