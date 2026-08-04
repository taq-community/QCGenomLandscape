test_that("parse_latlon handles cardinal directions", {
  expect_equal(parse_latlon("45.5 N 73.6 W"), list(lat = 45.5, lon = -73.6))
  expect_equal(parse_latlon("45.5N 73.6W"), list(lat = 45.5, lon = -73.6))
  expect_equal(parse_latlon("33.9 S 18.4 E"), list(lat = -33.9, lon = 18.4))
})

test_that("parse_latlon handles missing/malformed input", {
  expect_equal(parse_latlon(NA), list(lat = NA_real_, lon = NA_real_))
  expect_equal(parse_latlon(""), list(lat = NA_real_, lon = NA_real_))
  expect_equal(parse_latlon("45.5 N"), list(lat = NA_real_, lon = NA_real_))
})
