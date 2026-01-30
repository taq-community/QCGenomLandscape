bold_data <- read.delim(
  file = "results/bold_qc_data.tsv",
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  comment.char = ""
)

# Extract latitude and longitude from coord column
bold_data <- bold_data |>
  dplyr::mutate(
    latitude = as.numeric(stringr::str_extract(coord, "(?<=\\[)[0-9.-]+")),
    longitude = as.numeric(stringr::str_extract(coord, "(?<=, )[0-9.-]+(?=\\])"))
  )

bold_sf <- bold_data |>
    dplyr::filter(!is.na(latitude) & !is.na(longitude)) |>
    sf::st_as_sf(
        coords = c("longitude", "latitude"),
        crs = 4326
    )

mapview::mapview(bold_sf)

