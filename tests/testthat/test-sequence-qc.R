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

test_that("score_sequence_quality computes length/N%/GC% correctly", {
  result <- score_sequence_quality(c("acgtacgt", "acgtnnnn"))

  expect_equal(result$seq_length, c(8, 8))
  expect_equal(result$n_count, c(0, 4))
  expect_equal(result$n_pct, c(0, 50))
  expect_equal(result$gc_pct, c(50, 25))
})
