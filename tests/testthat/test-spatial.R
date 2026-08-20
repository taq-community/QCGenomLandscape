square_boundary <- function(xmax = 100000, ymax = 100000, crs = 3347) {
  sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(xmax, 0), c(xmax, ymax), c(0, ymax), c(0, 0)))),
    crs = crs
  ) |> sf::st_as_sf()
}

test_that("make_hex_grid counts points per cell and preserves the total", {
  boundary_sf <- square_boundary()
  points_sf <- sf::st_as_sf(
    data.frame(x = c(10000, 50000, 90000), y = c(10000, 50000, 90000)),
    coords = c("x", "y"), crs = 3347
  )

  hex <- make_hex_grid(points_sf, boundary_sf, cellsize = 20000, crs_proj = 3347, log_scale = TRUE)

  expect_s3_class(hex, "sf")
  expect_true(all(hex$n_records > 0))
  expect_equal(sum(hex$n_records), 3)
  expect_true("log_count" %in% names(hex))
})

test_that("make_hex_grid omits log_count when log_scale = FALSE", {
  boundary_sf <- square_boundary()
  points_sf <- sf::st_as_sf(data.frame(x = 50000, y = 50000), coords = c("x", "y"), crs = 3347)

  hex <- make_hex_grid(points_sf, boundary_sf, cellsize = 20000, crs_proj = 3347, log_scale = FALSE)

  expect_false("log_count" %in% names(hex))
})

test_that("flag_within_boundary flags points inside vs outside", {
  boundary_sf <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)))),
    crs = 4326
  ) |> sf::st_as_sf()

  points_sf <- sf::st_as_sf(
    data.frame(id = 1:2, x = c(5, 50), y = c(5, 50)),
    coords = c("x", "y"), crs = 4326
  )

  flagged <- flag_within_boundary(points_sf, boundary_sf)

  expect_equal(flagged$in_boundary, c(TRUE, FALSE))
})

test_that("flag_within_boundary reprojects the boundary when CRS differs", {
  boundary_sf <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)))),
    crs = 4326
  ) |> sf::st_as_sf() |> sf::st_transform(3857)

  points_sf <- sf::st_as_sf(data.frame(x = 5, y = 5), coords = c("x", "y"), crs = 4326)

  flagged <- flag_within_boundary(points_sf, boundary_sf, flag_name = "inside")

  expect_true(flagged$inside)
})

# ---- extract_edna_occurrences ----

test_that("extract_edna_occurrences normalizes NCBI + BOLD into a common shape", {
  ncbi_results <- tibble::tibble(
    organism = c("Alces alces", "Ursus americanus", "Missing coords"),
    lat_lon = c("45.5 N 73.6 W", "48.0 N 71.0 W", NA),
    collection_date = c("2020-06-15", "2019-01-01", NA),
    accession = c("ACC1", "ACC2", "ACC3")
  )

  bold_raw <- paste(
    "species\tcoord\tcollection_date_start\trecord_id",
    "Salmo salar\t[48.5, -68.5]\t2021-07-01\tBOLD1",
    "\t[45.0, -73.0]\t2021-01-01\tBOLD2", # no species -- unidentified, dropped
    "Esox lucius\t\t2021-01-01\tBOLD3", # no coord, dropped
    sep = "\n"
  )

  out <- extract_edna_occurrences(ncbi_results, bold_raw)

  expect_named(out, c("source", "species", "lon", "lat", "date", "record_id"))
  expect_setequal(out$record_id, c("ACC1", "ACC2", "BOLD1"))

  ncbi_row <- out[out$record_id == "ACC1", ]
  expect_equal(ncbi_row$source, "NCBI")
  expect_equal(ncbi_row$lat, 45.5)
  expect_equal(ncbi_row$lon, -73.6)
  expect_equal(ncbi_row$date, as.Date("2020-06-15"))

  bold_row <- out[out$record_id == "BOLD1", ]
  expect_equal(bold_row$source, "BOLD")
  expect_equal(bold_row$lat, 48.5)
  expect_equal(bold_row$lon, -68.5)
})
