species_df <- tibble::tibble(
  species = c("Alces alces", "Ursus americanus"),
  group_en = c("Mammals", "Mammals")
)
query_primers <- tibble::tibble(
  group = c("mammals"),
  query_marker = c("COI[Gene]")
)

test_that("build_ncbi_queries appends the voucher filter when voucher = TRUE", {
  queries <- build_ncbi_queries(species_df, query_primers, voucher = TRUE)

  expect_length(queries, 2)
  expect_true(all(grepl("AND voucher\\[Title\\]$", queries)))
  expect_true(any(grepl("^Alces alces\\[Organism\\]", queries)))
})

test_that("build_ncbi_queries omits the voucher filter when voucher = FALSE", {
  queries <- build_ncbi_queries(species_df, query_primers, voucher = FALSE)

  expect_length(queries, 2)
  expect_false(any(grepl("voucher", queries)))
})

test_that("build_ncbi_queries drops species with no matching primer marker", {
  species_no_match <- tibble::tibble(species = "Mystery species", group_en = "Unmapped Group")
  queries <- build_ncbi_queries(species_no_match, query_primers)

  expect_length(queries, 0)
})

# ---- fetch_ncbi_sequences: stubbed search_fn/summary_fn, no network ----

fake_search_fn <- function(db, term, retmax) {
  if (retmax == 0) {
    return(list(count = 2))
  }
  list(count = 2, ids = c("111", "222"))
}

fake_summary_fn <- function(db, id) {
  purrr::map(id, function(x) {
    list(
      uid = x,
      accessionversion = paste0("ACC", x),
      title = "Test record",
      taxid = "9999",
      organism = "Alces alces",
      moltype = "genomic DNA",
      topology = "linear",
      genome = "genomic",
      slen = "658",
      createdate = "2020/01/01",
      updatedate = "2020/01/01",
      subtype = "specimen_voucher|country",
      subname = "ABC:123|Canada: Quebec"
    )
  })
}

test_that("fetch_ncbi_sequences builds a results tibble from stubbed search/summary functions", {
  out <- fetch_ncbi_sequences(
    queries = "Alces alces[Organism] AND COI[Gene]",
    search_fn = fake_search_fn,
    summary_fn = fake_summary_fn,
    progress = FALSE
  )

  expect_named(out, c("results", "deficient_queries", "high_id_queries"))
  expect_equal(nrow(out$results), 2)
  expect_equal(out$results$accession, c("ACC111", "ACC222"))
  expect_equal(out$results$specimen_voucher, c("ABC:123", "ABC:123"))
  expect_length(out$deficient_queries, 0)
  expect_length(out$high_id_queries, 0)
})

test_that("fetch_ncbi_sequences flags queries above the high-ID threshold", {
  out <- fetch_ncbi_sequences(
    queries = "Alces alces[Organism] AND COI[Gene]",
    search_fn = fake_search_fn,
    summary_fn = fake_summary_fn,
    high_id_threshold = 1,
    progress = FALSE
  )

  expect_length(out$high_id_queries, 1)
  expect_equal(out$high_id_queries[[1]]$id_count, 2)
})

test_that("fetch_ncbi_sequences records a deficient query when search_fn errors", {
  erroring_search_fn <- function(db, term, retmax) stop("simulated NCBI error")

  out <- fetch_ncbi_sequences(
    queries = "Alces alces[Organism] AND COI[Gene]",
    search_fn = erroring_search_fn,
    summary_fn = fake_summary_fn,
    progress = FALSE
  )

  expect_equal(nrow(out$results), 0)
  expect_length(out$deficient_queries, 1)
  expect_equal(out$deficient_queries[[1]]$error_type, "entrez_search_count")
})

test_that("fetch_ncbi_sequences returns an empty result for a zero-count query", {
  zero_count_search_fn <- function(db, term, retmax) list(count = 0)

  out <- fetch_ncbi_sequences(
    queries = "Nonexistent species[Organism]",
    search_fn = zero_count_search_fn,
    summary_fn = fake_summary_fn,
    progress = FALSE
  )

  expect_equal(nrow(out$results), 0)
  expect_length(out$deficient_queries, 0)
})
