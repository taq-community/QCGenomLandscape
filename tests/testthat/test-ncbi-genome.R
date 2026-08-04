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

test_that("fetch_ncbi_genomes takes the fast zero-hit path without per-species calls", {
  call_count <- 0
  fake_search_fn <- function(db, term, retmax) {
    call_count <<- call_count + 1
    list(count = 0, ids = character(0))
  }

  result <- fetch_ncbi_genomes(
    c("Species A", "Species B", "Species C"),
    batch_size = 25, progress = FALSE, search_fn = fake_search_fn
  )

  expect_equal(nrow(result), 3)
  expect_false(any(result$nuclear_genome))
  expect_false(any(result$mitochondrial_genome))
  # exactly 2 calls total: one genome existence check, one mito existence
  # check -- NOT 6 (2 per species) as the unbatched version would make
  expect_equal(call_count, 2)
})

test_that("fetch_ncbi_genomes falls back to per-species nuclear resolution when the batch has a hit", {
  fake_search_fn <- function(db, term, retmax) {
    is_batched <- grepl(" OR ", term, fixed = TRUE)

    if (is_batched) {
      # batch existence check: genome db has a hit, mito db doesn't
      if (db == "genome") return(list(count = 1, ids = "SOME_ID"))
      return(list(count = 0, ids = character(0)))
    }

    # per-species fallback query -- only "Species A" has a genome
    if (db == "genome" && grepl("Species A\\[Organism\\]", term)) {
      return(list(count = 1, ids = "GENOME_A"))
    }
    list(count = 0, ids = character(0))
  }

  result <- fetch_ncbi_genomes(
    c("Species A", "Species B"),
    batch_size = 25, progress = FALSE, search_fn = fake_search_fn
  )

  expect_equal(nrow(result), 2)
  expect_true(result$nuclear_genome[result$species == "Species A"])
  expect_false(result$nuclear_genome[result$species == "Species B"])
})

test_that("fetch_ncbi_genomes attributes a batched mitochondrial hit back to the right species", {
  fake_search_fn <- function(db, term, retmax) {
    if (db == "genome") return(list(count = 0, ids = character(0))) # no nuclear hits
    list(count = 2, ids = c("M1", "M2"))
  }
  fake_summary_fn <- function(db, id) {
    purrr::map(id, function(x) {
      organism <- if (x == "M1") "Species A" else "Species C" # M2 doesn't match any batch species
      list(uid = x, organism = organism)
    })
  }

  result <- fetch_ncbi_genomes(
    c("Species A", "Species B"),
    batch_size = 25, progress = FALSE, search_fn = fake_search_fn, summary_fn = fake_summary_fn
  )

  expect_equal(nrow(result), 2)
  expect_true(result$mitochondrial_genome[result$species == "Species A"])
  expect_equal(result$mitochondrial_accessions[result$species == "Species A"], "M1")
  expect_false(result$mitochondrial_genome[result$species == "Species B"])
})

test_that("fetch_ncbi_genomes never does per-species mitochondrial calls, even with a hit", {
  # the whole point of the fix: mito resolution is always exactly 1 search
  # call per batch (no existence pre-check, no per-species fallback),
  # regardless of whether it hits
  mito_call_count <- 0
  fake_search_fn <- function(db, term, retmax) {
    if (db == "nucleotide") mito_call_count <<- mito_call_count + 1
    if (db == "genome") return(list(count = 0, ids = character(0)))
    list(count = 3, ids = c("M1", "M2", "M3"))
  }
  fake_summary_fn <- function(db, id) {
    purrr::map(id, function(x) list(uid = x, organism = "Species A"))
  }

  fetch_ncbi_genomes(
    c("Species A", "Species B", "Species C"),
    batch_size = 25, progress = FALSE, search_fn = fake_search_fn, summary_fn = fake_summary_fn
  )

  expect_equal(mito_call_count, 1)
})

test_that("fetch_ncbi_genomes keeps nuclear and mitochondrial resolution independent", {
  # regression test for the actual bug: a mito-only hit must not drag the
  # whole batch into per-species nuclear resolution
  nuclear_calls <- 0
  fake_search_fn <- function(db, term, retmax) {
    if (db == "genome") {
      nuclear_calls <<- nuclear_calls + 1
      return(list(count = 0, ids = character(0))) # never a nuclear hit
    }
    list(count = 1, ids = "M1") # always a mito hit
  }
  fake_summary_fn <- function(db, id) list(list(uid = "M1", organism = "Species A"))

  fetch_ncbi_genomes(
    c("Species A", "Species B", "Species C"),
    batch_size = 25, progress = FALSE, search_fn = fake_search_fn, summary_fn = fake_summary_fn
  )

  # exactly 1 nuclear existence check for the whole batch -- NOT 3
  # (1 per species), which is what the pre-fix code did whenever mito hit
  expect_equal(nuclear_calls, 1)
})

test_that("fetch_ncbi_genomes splits species into multiple batches", {
  genome_terms <- character(0)
  fake_search_fn <- function(db, term, retmax) {
    # the mito query's own fixed clauses always contain " OR " regardless of
    # batch size, so only `genome`-db terms are a clean multi-species signal
    if (db == "genome") genome_terms <<- c(genome_terms, term)
    list(count = 0, ids = character(0))
  }

  fetch_ncbi_genomes(
    c("Species A", "Species B", "Species C"),
    batch_size = 2, progress = FALSE, search_fn = fake_search_fn
  )

  # batch 1: "(Species A[Organism] OR Species B[Organism])"; batch 2: just C
  expect_equal(length(genome_terms), 2)
  expect_true(grepl(" OR ", genome_terms[1], fixed = TRUE))
  expect_false(grepl(" OR ", genome_terms[2], fixed = TRUE))
})

test_that("fetch_ncbi_genomes deduplicates species", {
  fake_search_fn <- function(db, term, retmax) list(count = 0, ids = character(0))

  result <- fetch_ncbi_genomes(
    c("Species A", "Species A", "Species B"),
    batch_size = 25, progress = FALSE, search_fn = fake_search_fn
  )

  expect_equal(nrow(result), 2)
})
