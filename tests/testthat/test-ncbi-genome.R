test_that("query_full_genome reports nuclear and mitochondrial genome availability", {
  fake_search_fn <- function(db, term, retmax) {
    if (db == "genome") {
      return(list(count = 1, ids = "GENOME123"))
    }
    list(count = 2, ids = c("MITO1", "MITO2"))
  }

  result <- query_full_genome("Alces alces", search_fn = fake_search_fn)

  expect_true(result$nuclear_genome)
  expect_true(result$mitochondrial_genome)
  expect_equal(result$nuclear_count, 1)
  expect_equal(result$mitochondrial_count, 2)
  expect_equal(result$nuclear_accessions, "GENOME123")
  expect_equal(result$mitochondrial_accessions, "MITO1,MITO2")
})

test_that("query_full_genome reports no genome when counts are zero", {
  fake_search_fn <- function(db, term, retmax) list(count = 0, ids = character(0))

  result <- query_full_genome("Nonexistent species", search_fn = fake_search_fn)

  expect_false(result$nuclear_genome)
  expect_false(result$mitochondrial_genome)
})

test_that("query_full_genome captures errors without throwing", {
  erroring_search_fn <- function(db, term, retmax) stop("simulated NCBI error")

  result <- query_full_genome("Alces alces", search_fn = erroring_search_fn)

  expect_equal(result$error, "simulated NCBI error")
})

test_that("fetch_ncbi_genomes maps query_full_genome over a species vector", {
  fake_search_fn <- function(db, term, retmax) list(count = 0, ids = character(0))

  result <- fetch_ncbi_genomes(c("Species A", "Species B"), progress = FALSE, search_fn = fake_search_fn)

  expect_equal(nrow(result), 2)
  expect_equal(result$species, c("Species A", "Species B"))
})
