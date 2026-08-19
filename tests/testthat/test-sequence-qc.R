test_that("has_stop_codon_coi detects an in-frame stop codon", {
  # frame 1 codons: aaa | taa | aaa -- "taa" is a stop in the invertebrate
  # mitochondrial code (SGC4)
  expect_true(has_stop_codon_coi("aaataaaaa"))
})

test_that("has_stop_codon_coi returns FALSE when no frame has a stop codon", {
  # homopolymeric sequence: every codon in every frame is "aaa" (Lys), never a stop
  expect_false(has_stop_codon_coi("aaaaaaaaa"))
})

test_that("has_stop_codon_coi returns NA for NA input", {
  expect_true(is.na(has_stop_codon_coi(NA)))
})

test_that("has_stop_codon_coi is vectorized and preserves position, including NA and short sequences", {
  result <- has_stop_codon_coi(c("aaataaaaa", "aaaaaaaaa", NA, "aa"))

  expect_equal(result, c(TRUE, FALSE, NA, FALSE))
})

test_that("has_stop_codon_coi gives the same answer vectorized as one at a time", {
  seqs <- c("aaataaaaa", "aaaaaaaaa", "atgtaaatg", "acgtacgtacgt")

  vectorized <- has_stop_codon_coi(seqs)
  scalar <- vapply(seqs, has_stop_codon_coi, logical(1), USE.NAMES = FALSE)

  expect_equal(vectorized, scalar)
})

test_that("score_sequence_quality computes length/N%/GC% correctly", {
  result <- score_sequence_quality(c("acgtacgt", "acgtnnnn"))

  expect_equal(result$seq_length, c(8, 8))
  expect_equal(result$n_count, c(0, 4))
  expect_equal(result$n_pct, c(0, 50))
  expect_equal(result$gc_pct, c(50, 25))
})

test_that("build_sequence_qc_table parses batches, scores quality, and flags stop codons", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )
  # two "batches", mirroring fetch_gb_records()'s list-of-blobs shape
  gb_records <- list(gb_text, gb_text)

  out <- build_sequence_qc_table(
    gb_records,
    log_every_batch = 1, log_every_seq = 1, progress = FALSE, cores = 2
  )

  expect_equal(nrow(out), 6)
  expect_named(
    out,
    c("accession", "gene", "sequence", "seq_length", "n_count", "n_pct", "gc_pct", "has_stop")
  )
  expect_true(is.na(out$sequence[out$accession == "GHI789"][1]))
  expect_true(is.na(out$has_stop[out$accession == "GHI789"][1]))
  expect_false(is.na(out$has_stop[out$accession == "ABC123"][1]))
})

test_that("build_sequence_qc_table drops a batch that fails to parse instead of erroring", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )
  # NULL makes parse_gb_records() error (non-character argument) --
  # mclapply() should catch it, not propagate it
  gb_records <- list(gb_text, NULL)

  out <- build_sequence_qc_table(
    gb_records,
    log_every_batch = 1, log_every_seq = 1, progress = FALSE, cores = 2
  )

  expect_equal(nrow(out), 3)
  expect_setequal(out$accession, c("ABC123", "DEF456", "GHI789"))
})
