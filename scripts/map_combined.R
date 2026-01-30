library(shiny)
library(bslib)
library(leaflet)
library(sf)
library(dplyr)
library(stringr)
library(purrr)

# -----------------------------------------------------------------------------
# Load utility functions
# -----------------------------------------------------------------------------
source("scripts/utils.R"))

# -----------------------------------------------------------------------------
# Data Loading and Preprocessing
# -----------------------------------------------------------------------------

# Load taxonomic reference data
bb_taxo <- read.csv("data/bdqc_list_01122025.csv") |>
  select(valid_scientific_name, rank, kingdom, phylum, class, order, family, genus, species) |>
  filter(rank == "species") |>
  distinct()

# Load BOLD data
bold_data <- read.delim(
    file = "results/bold_qc_data.tsv",
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    comment.char = ""
) |>
    filter(identification_rank == "species") |>
    mutate(
        latitude = as.numeric(str_extract(coord, "(?<=\\[)[0-9.-]+")),
        longitude = as.numeric(str_extract(coord, "(?<=, )[0-9.-]+(?=\\])")),
        date = as.Date(collection_date_start),
        source = "BOLD",
        species = species
    ) |>
    select(
        species,
        latitude, longitude, date, source
    ) |>
    filter(!is.na(latitude) & !is.na(longitude)) |>
    left_join(
        bb_taxo,
        by = "species"
    ) |>
    filter(!is.na(valid_scientific_name))


# Load NCBI data
ncbi_data <- readRDS("results/ncbi_results.rds") |>
  mutate(
    parsed_coords = map(lat_lon, parse_latlon),
    latitude = map_dbl(parsed_coords, "lat"),
    longitude = map_dbl(parsed_coords, "lon"),
    date = as.Date(collection_date, format = "%d-%b-%Y"),
    source = "NCBI",
    identification_rank = NA_character_,
    species = str_extract(query, "^\\w+\\s+\\w+")
  ) |>
  select(-parsed_coords) |>
  left_join(
    bb_taxo,
    by = "species"
  ) |>
  select(
    valid_scientific_name, rank, kingdom, phylum, class, order, family, genus, species,
    identification_rank, latitude, longitude, date, source
  ) |>
  filter(!is.na(latitude) & !is.na(longitude))

# Combine datasets
combined_data <- bind_rows(bold_data, ncbi_data)

# Load Quebec boundary
qc <- read_sf("data/canvec_1M_CA_Admin.gdb", layer = "geo_political_region_2") |>
  filter(jurisdiction == 102) |>
  st_transform(4326)

# Filter to Quebec
combined_sf <- st_as_sf(
  combined_data,
  coords = c("longitude", "latitude"),
  crs = 4326
)
combined_sf$in_qc <- lengths(st_within(combined_sf$geometry, qc)) > 0
combined_qc <- combined_sf |> filter(in_qc)

# Extract coordinates back for filtering
combined_qc <- combined_qc |>
  mutate(
    longitude = st_coordinates(geometry)[, 1],
    latitude = st_coordinates(geometry)[, 2]
  )

mapview::mapview(combined_qc)

# Get unique values for filters
families <- sort(unique(combined_qc$family[!is.na(combined_qc$family)]))
ranks <- sort(unique(combined_qc$identification_rank[!is.na(combined_qc$identification_rank)]))
date_range <- range(combined_qc$date, na.rm = TRUE)

# Quebec projected CRS for hexagon creation
qc_proj <- st_transform(qc, 3347)

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

ui <- page_sidebar(
  title = "Quebec Genomic Landscape - NCBI & BOLD Data",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 300,
    selectInput(
      "family",
      "Family",
      choices = c("All" = "", families),
      selected = "",
      multiple = TRUE
    ),
    selectInput(
      "rank",
      "Identification Rank",
      choices = c("All" = "", ranks),
      selected = "",
      multiple = TRUE
    ),
    sliderInput(
      "date_range",
      "Collection Date Range",
      min = date_range[1],
      max = date_range[2],
      value = date_range,
      timeFormat = "%Y-%m-%d"
    ),
    checkboxInput(
      "use_log",
      "Use log(count) for hexagon coloring",
      value = FALSE
    ),
    hr(),
    p(
      class = "text-muted",
      paste0(
        "Total records: ", nrow(combined_qc),
        " (BOLD: ", sum(combined_qc$source == "BOLD"),
        ", NCBI: ", sum(combined_qc$source == "NCBI"), ")"
      )
    )
  ),
  card(
    card_header("Occurrence Map"),
    leafletOutput("map", height = "calc(100vh - 200px)")
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {
  # Reactive filtered data

filtered_data <- reactive({
    data <- combined_qc

    # Filter by family
    if (length(input$family) > 0 && !all(input$family == "")) {
      data <- data |> filter(family %in% input$family)
    }

    # Filter by rank
    if (length(input$rank) > 0 && !all(input$rank == "")) {
      data <- data |> filter(identification_rank %in% input$rank)
    }

    # Filter by date
    if (!is.null(input$date_range)) {
      data <- data |>
        filter(
          is.na(date) |
            (date >= input$date_range[1] & date <= input$date_range[2])
        )
    }

    data
  })

  # Create hexagon grid with counts
  hex_data <- reactive({
    data <- filtered_data()

    if (nrow(data) == 0) {
      return(NULL)
    }

    # Project data
    data_proj <- st_transform(data, 3347)

    # Create hexagonal grid over Quebec (10km = 10000m cell size)
    hex_grid <- st_make_grid(
      qc_proj,
      cellsize = 10000,
      square = FALSE
    ) |>
      st_as_sf() |>
      st_intersection(qc_proj)

    # Count points in each hexagon
    hex_grid$n_records <- lengths(st_intersects(hex_grid, data_proj))

    # Filter to hexagons with at least one record
    hex_with_data <- hex_grid |> filter(n_records > 0)

    if (nrow(hex_with_data) == 0) {
      return(NULL)
    }

    # Add log count
    hex_with_data <- hex_with_data |>
      mutate(log_count = log(n_records + 1))

    # Transform back to WGS84 for mapping
    st_transform(hex_with_data, 4326)
  })

  # Render map
  output$map <- renderLeaflet({
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(lng = -72, lat = 52, zoom = 5)
  })

  # Update hexagons when data changes
  observe({
    hex <- hex_data()

    leafletProxy("map") |>
      clearShapes() |>
      clearControls()

    if (is.null(hex) || nrow(hex) == 0) {
      return()
    }

    # Choose color variable based on checkbox
    if (input$use_log) {
      color_var <- hex$log_count
      legend_title <- "log(Count)"
    } else {
      color_var <- hex$n_records
      legend_title <- "Count"
    }

    # Create color palette
    pal <- colorNumeric(
      palette = "YlOrRd",
      domain = color_var
    )

    leafletProxy("map") |>
      addPolygons(
        data = hex,
        fillColor = ~pal(color_var),
        fillOpacity = 0.7,
        color = "#444444",
        weight = 0.5,
        popup = ~paste0(
          "<strong>Records: </strong>", n_records, "<br>",
          "<strong>log(Count): </strong>", round(log_count, 2)
        )
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = color_var,
        title = legend_title,
        opacity = 0.7
      )
  })
}

# -----------------------------------------------------------------------------
# Run App
# -----------------------------------------------------------------------------

shinyApp(ui, server)
