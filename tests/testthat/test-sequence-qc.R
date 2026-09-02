test_that("has_stop_codon_coi returns FALSE (not flagged) when at least one frame is clean", {
  # frame 1 codons: aaa | taa | aaa -- "taa" is a stop, but frames 2/3
  # ("aataaaaa"/"ataaaaa" shifted) are stop-free -- one clean frame is
  # enough to not flag the sequence, since that's presumably the real one
  expect_false(has_stop_codon_coi("aaataaaaa"))
  # homopolymer: every codon in every frame is "aaa" (Lys) -- all 3 clean
  expect_false(has_stop_codon_coi("aaaaaaaaa"))
})

test_that("has_stop_codon_coi returns TRUE when no frame is clean", {
  # verified (scripts/... brute-force search): all 3 reading frames of this
  # 30bp sequence hit a stop codon somewhere -- the rare, meaningful signal
  expect_true(has_stop_codon_coi("tgaggcgtatagaaacttagcgccgtagca"))
})

test_that("has_stop_codon_coi returns NA for NA input", {
  expect_true(is.na(has_stop_codon_coi(NA)))
})

test_that("has_stop_codon_coi is vectorized and preserves position, including NA and short sequences", {
  result <- has_stop_codon_coi(c(
    "aaataaaaa", "tgaggcgtatagaaacttagcgccgtagca", NA, "aa"
  ))

  expect_equal(result, c(FALSE, TRUE, NA, FALSE))
})

test_that("has_stop_codon_coi gives the same answer vectorized as one at a time", {
  seqs <- c("aaataaaaa", "aaaaaaaaa", "atgtaaatg", "acgtacgtacgt", "tgaggcgtatagaaacttagcgccgtagca")

  vectorized <- has_stop_codon_coi(seqs)
  scalar <- vapply(seqs, has_stop_codon_coi, logical(1), USE.NAMES = FALSE)

  expect_equal(vectorized, scalar)
})

test_that("flag_length_outliers flags a length far from its group's median", {
  lengths <- c(640, 655, 648, 651, 645, 652, 649, 647, 653, 646, 60)
  group <- rep("COI", length(lengths))

  out <- flag_length_outliers(lengths, group, min_n = 5)

  expect_equal(out, c(rep(FALSE, 10), TRUE))
})

test_that("flag_length_outliers evaluates each group independently", {
  lengths <- c(640, 655, 648, 651, 645, 1400, 1420, 1390, 1410, 1405)
  group <- rep(c("COI", "rbcL"), each = 5)

  out <- flag_length_outliers(lengths, group, min_n = 5)

  expect_true(all(!out))
})

test_that("flag_length_outliers leaves groups smaller than min_n as NA", {
  out <- flag_length_outliers(c(640, 5), c("COI", "COI"), min_n = 10)

  expect_true(all(is.na(out)))
})

test_that("flag_length_outliers falls back to exact-match when a group has zero variance", {
  out <- flag_length_outliers(c(650, 650, 650, 650, 650, 620), rep("COI", 6), min_n = 5)

  expect_equal(out, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
})

test_that("flag_length_outliers propagates NA seq_length without breaking the group", {
  out <- flag_length_outliers(c(640, 655, 648, 651, 645, NA), rep("COI", 6), min_n = 5)

  expect_true(is.na(out[6]))
  expect_false(any(out[1:5]))
})

test_that("flag_barcode_gap_outliers flags a sequence far from its group's usual distance", {
  seqs <- c(
    "ATGAAATTTGGGCCCAAATTTGGGCCC",
    "ATGAAATTTGGGCCCAAATTTGGGCCT",
    "ATGAAATTTGGGCCCAAATTTGGGCCC",
    "TTTTTTTTTTTTTTTTTTTTTTTTTTT"
  )
  species <- rep("Alces alces", 4)
  group <- rep("COI", 4)

  out <- flag_barcode_gap_outliers(seqs, species, group, min_n = 3, cores = 2)

  expect_equal(out, c(FALSE, FALSE, FALSE, TRUE))
})

test_that("flag_barcode_gap_outliers only compares within the same (species, group)", {
  seqs <- c("ATGAAATTTGGGCCC", "ATGAAATTTGGGCCT", "TTTTTTTTTTTTTTT")
  species <- c("Alces alces", "Alces alces", "Ursus americanus")
  group <- c("COI", "COI", "COI")

  # Ursus's lone sequence is in a group of size 1 -- below min_n, left NA,
  # not compared against Alces's sequences despite being very different
  out <- flag_barcode_gap_outliers(seqs, species, group, min_n = 2, cores = 2)

  expect_equal(out[1:2], c(FALSE, FALSE))
  expect_true(is.na(out[3]))
})

test_that("flag_barcode_gap_outliers leaves groups outside [min_n, max_n] as NA", {
  seqs <- c("ATGAAATTTGGGCCC", "ATGAAATTTGGGCCT")
  out <- flag_barcode_gap_outliers(seqs, c("A", "A"), c("COI", "COI"), min_n = 5, cores = 2)

  expect_true(all(is.na(out)))
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
    log_every_batch = 1, progress = FALSE, cores = 2
  )

  expect_equal(nrow(out), 6)
  expect_named(
    out,
    c(
      "accession", "definition", "gene", "sequence", "seq_length", "n_count", "n_pct", "gc_pct",
      "is_complete_genome", "has_stop"
    )
  )
  expect_true(is.na(out$sequence[out$accession == "GHI789"][1]))
  expect_true(is.na(out$has_stop[out$accession == "GHI789"][1]))
  expect_false(is.na(out$has_stop[out$accession == "ABC123"][1]))
  # none of the fixture's DEFINITION lines say "complete genome"
  expect_false(any(out$is_complete_genome))
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
    log_every_batch = 1, progress = FALSE, cores = 2
  )

  expect_equal(nrow(out), 3)
  expect_setequal(out$accession, c("ABC123", "DEF456", "GHI789"))
})
