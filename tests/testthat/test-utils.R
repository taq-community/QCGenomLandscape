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

test_that("parse_bold_coord parses '[lat, lon]' strings, including negative longitudes", {
  expect_equal(parse_bold_coord("[45.45, -75.77]"), list(lat = 45.45, lon = -75.77))
  expect_equal(parse_bold_coord("[46.25, -72]"), list(lat = 46.25, lon = -72))
})

test_that("parse_bold_coord handles missing/malformed input", {
  expect_equal(parse_bold_coord(NA), list(lat = NA_real_, lon = NA_real_))
  expect_equal(parse_bold_coord(""), list(lat = NA_real_, lon = NA_real_))
  expect_equal(parse_bold_coord("[45.45]"), list(lat = NA_real_, lon = NA_real_))
})

test_that("parse_gb_collection_date handles every known GenBank/BOLD date precision", {
  out <- parse_gb_collection_date(c(
    "28-Sep-2014", "Sep-2023", "2012-08", "2020", "2020-06-15"
  ))

  expect_equal(out, as.Date(c(
    "2014-09-28", "2023-09-01", "2012-08-01", "2020-01-01", "2020-06-15"
  )))
})

test_that("parse_gb_collection_date returns NA for junk, ranges, and missing values", {
  out <- parse_gb_collection_date(c(
    "vicinity of Seventy Five Mile Creek", "58.298 N 8.539 E",
    "2012/2014", "May-2012/May-2017", NA, ""
  ))

  expect_true(all(is.na(out)))
})

test_that("parse_gb_collection_date never errors on a vector with no matches at all", {
  # as.Date()'s default tryFormats throws in exactly this situation --
  # regression guard for that base-R quirk
  expect_no_error(parse_gb_collection_date(c("junk1", "junk2", "junk3")))
})
