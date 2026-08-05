species_df <- tibble::tibble(
  species = c("Alces alces", "Ursus americanus"),
  group_en = c("Mammals", "Mammals")
)
query_primers <- tibble::tibble(
  group = c("mammals"),
  query_marker = c("COI[Gene]")
)

test_that("build_ncbi_queries builds one query per species/marker pair", {
  queries <- build_ncbi_queries(species_df, query_primers)

  expect_length(queries, 2)
  expect_false(any(grepl("voucher", queries)))
  expect_true(any(grepl("^Alces alces\\[Organism\\]", queries)))
})

test_that("build_ncbi_queries drops species with no matching primer marker", {
  species_no_match <- tibble::tibble(species = "Mystery species", group_en = "Unmapped Group")
  queries <- build_ncbi_queries(species_no_match, query_primers)

  expect_length(queries, 0)
})

test_that("build_ncbi_queries OR-combines species sharing a marker when batch_size > 1", {
  queries <- build_ncbi_queries(species_df, query_primers, batch_size = 25)

  expect_length(queries, 1)
  expect_true(grepl("^\\(Alces alces\\[Organism\\] OR Ursus americanus\\[Organism\\]\\)", queries))
})

test_that("build_ncbi_queries splits into multiple batches once batch_size is exceeded", {
  queries <- build_ncbi_queries(species_df, query_primers, batch_size = 1)

  expect_length(queries, 2)
  expect_false(any(grepl("OR", queries)))
})

test_that("build_ncbi_queries never mixes species across different markers in one batch", {
  species_two_groups <- tibble::tibble(
    species = c("Alces alces", "Salmo salar"),
    group_en = c("Mammals", "Fish")
  )
  primers_two_groups <- tibble::tibble(
    group = c("mammals", "fish"),
    query_marker = c("COI[Gene]", "cytb[Gene]")
  )

  queries <- build_ncbi_queries(species_two_groups, primers_two_groups, batch_size = 25)

  expect_length(queries, 2)
  expect_false(any(grepl("OR", queries)))
})

# ---- fetch_ncbi_sequences: stubbed search_fn/summary_fn, no network ----

fake_search_fn <- function(db, term, retmax) {
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
  expect_equal(out$results$is_voucher, c(FALSE, FALSE))
  expect_length(out$deficient_queries, 0)
  expect_length(out$high_id_queries, 0)
})

test_that("fetch_ncbi_sequences infers is_voucher from the record title, case-insensitively", {
  mixed_title_summary_fn <- function(db, id) {
    titles <- c(
      "111" = "Alces alces voucher specimen COI gene, partial cds",
      "222" = "Alces alces VOUCHER XYZ COI gene",
      "333" = "Alces alces isolate ABC COI gene, partial cds"
    )
    purrr::map(id, function(x) {
      list(
        uid = x,
        accessionversion = paste0("ACC", x),
        title = titles[[x]],
        organism = "Alces alces"
      )
    })
  }

  out <- fetch_ncbi_sequences(
    queries = "Alces alces[Organism] AND COI[Gene]",
    search_fn = function(db, term, retmax) list(count = 3, ids = c("111", "222", "333")),
    summary_fn = mixed_title_summary_fn,
    progress = FALSE
  )

  expect_equal(out$results$is_voucher, c(TRUE, TRUE, FALSE))
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
  expect_equal(out$deficient_queries[[1]]$error_type, "entrez_search_ids")
})

test_that("fetch_ncbi_sequences returns an empty result for a zero-count query", {
  zero_count_search_fn <- function(db, term, retmax) list(count = 0, ids = character(0))

  out <- fetch_ncbi_sequences(
    queries = "Nonexistent species[Organism]",
    search_fn = zero_count_search_fn,
    summary_fn = fake_summary_fn,
    progress = FALSE
  )

  expect_equal(nrow(out$results), 0)
  expect_length(out$deficient_queries, 0)
})
