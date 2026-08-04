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
