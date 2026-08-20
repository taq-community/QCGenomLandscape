test_that("atlas_molecular_datasets returns the expected fixed list", {
  ds <- atlas_molecular_datasets()

  expect_type(ds, "character")
  expect_true("International Barcode of Life project (iBOL)" %in% ds)
  expect_length(ds, 16)
})

test_that("build_hex_grid_cells returns every cell with a stable cell_id, not just occupied ones", {
  boundary_sf <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
    crs = 4326
  ) |> sf::st_as_sf()

  hex <- build_hex_grid_cells(boundary_sf, cellsize = 20000, crs_proj = 3857)

  expect_s3_class(hex, "sf")
  expect_true("cell_id" %in% names(hex))
  expect_equal(hex$cell_id, seq_len(nrow(hex)))
  expect_true(nrow(hex) > 1)
})

test_that("summarize_edna_contribution flags species-, cell-, and year-level novelty correctly", {
  # species A: eDNA and traditional share the same (cell, year) -> nothing novel
  # species B: eDNA only, no traditional record anywhere -> fully novel
  edna_agg <- tibble::tibble(
    species = c("Species A", "Species B"),
    cell_id = c(1L, 2L),
    year_obs = c(2020L, 2021L),
    n = c(1L, 1L)
  )
  traditional_agg <- tibble::tibble(
    species = "Species A",
    cell_id = 1L,
    year_obs = 2020L,
    n = 5L
  )

  out <- summarize_edna_contribution(edna_agg, traditional_agg)

  expect_equal(out$edna_only_species, "Species B")
  expect_equal(out$n_species_edna, 2)
  expect_equal(out$n_species_edna_only, 1)
  expect_equal(out$novel_cells$species, "Species B")
  expect_equal(out$novel_years$species, "Species B")
  expect_equal(out$pct_cells_novel, 50)
  expect_equal(out$pct_years_novel, 50)

  # species A is seen by both methods, so it's neither eDNA-only nor
  # traditional-only -- traditional has no species of its own here
  expect_equal(out$traditional_only_species, character(0))
  expect_equal(out$shared_species, "Species A")
  expect_equal(out$n_species_traditional, 1)
  expect_equal(out$n_species_traditional_only, 0)
  expect_equal(out$pct_species_traditional_only, 0)
  expect_equal(out$n_species_shared, 1)
  # shared = 1 (Species A), union = 2 (Species A, Species B)
  expect_equal(out$jaccard_species, 0.5)
})

test_that("fetch_atlas_parquet reuses an already-downloaded file instead of re-downloading", {
  dest_dir <- withr::local_tempdir()
  existing <- file.path(dest_dir, "atlas_public_2024-01-01.parquet")
  writeLines("placeholder", existing)

  called <- FALSE
  path <- fetch_atlas_parquet(
    export_date = "2024-01-01",
    dest_dir = dest_dir,
    download_fn = function(url, destfile) called <<- TRUE
  )

  expect_equal(path, existing)
  expect_false(called)
})

test_that("fetch_atlas_parquet downloads to dest_dir and picks the latest date when unset", {
  dest_dir <- withr::local_tempdir()
  downloaded <- list()

  path <- fetch_atlas_parquet(
    dest_dir = dest_dir,
    list_dates_fn = function(url) c("2024-01-01", "2025-06-15", "2024-12-01"),
    download_fn = function(url, destfile) {
      downloaded[[1]] <<- list(url = url, destfile = destfile)
      writeLines("placeholder", destfile)
    }
  )

  expect_equal(path, file.path(dest_dir, "atlas_public_2025-06-15.parquet"))
  expect_true(grepl("atlas_public_2025-06-15\\.parquet$", downloaded[[1]]$url))
})

test_that("compare_edna_atlas_coverage classifies molecular vs traditional and merges in our own eDNA pull", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("duckdb")

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;")

  # Species A: an iBOL (molecular) record and a traditional (eBird) record at
  # the SAME point/year -- should end up fully redundant, not novel.
  # Species B: an iBOL record with no traditional counterpart at all -- fully novel.
  DBI::dbWriteTable(con, "tmp_atlas", data.frame(
    valid_scientific_name = c("Species A", "Species A", "Species B"),
    longitude = c(0.5, 0.5, 0.9),
    latitude = c(0.5, 0.5, 0.9),
    year_obs = c(2020L, 2020L, 2021L),
    dataset_name = c(
      "International Barcode of Life project (iBOL)",
      "eBird Observation Dataset",
      "International Barcode of Life project (iBOL)"
    ),
    stringsAsFactors = FALSE
  ))
  parquet_path <- withr::local_tempfile(fileext = ".parquet")
  DBI::dbExecute(con, sprintf("COPY tmp_atlas TO '%s' (FORMAT PARQUET)", parquet_path))

  boundary_sf <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
    crs = 4326
  ) |> sf::st_as_sf()
  hex_grid <- build_hex_grid_cells(boundary_sf, cellsize = 20000, crs_proj = 3857)

  # Species C: only in our own extraction, nowhere in the Atlas fixture
  edna_occurrences <- tibble::tibble(
    species = "Species C",
    lon = 0.2,
    lat = 0.2,
    date = as.Date("2022-01-01")
  )

  out <- compare_edna_atlas_coverage(parquet_path, hex_grid, edna_occurrences, con = con)

  expect_setequal(unique(out$edna_agg$species), c("Species A", "Species B", "Species C"))
  expect_equal(unique(out$traditional_agg$species), "Species A")

  contribution <- summarize_edna_contribution(out$edna_agg, out$traditional_agg)
  expect_setequal(contribution$edna_only_species, c("Species B", "Species C"))
})
