test_that("parse_gb_records parses accession, gene, and sequence per record", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_equal(nrow(parsed), 3)
  expect_equal(parsed$accession, c("ABC123", "DEF456", "GHI789"))
})

test_that("parse_gb_records de-duplicates repeated /gene= tags within a record", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_equal(parsed$gene[parsed$accession == "ABC123"], "COI")
})

test_that("parse_gb_records collapses multiple distinct /gene= tags with ';'", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_equal(parsed$gene[parsed$accession == "DEF456"], "COI;Cytb")
})

test_that("parse_gb_records extracts only acgt characters from ORIGIN", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)
  seq1 <- parsed$sequence[parsed$accession == "ABC123"]

  expect_true(grepl("^[acgtACGT]+$", seq1))
  expect_equal(seq1, "acgtacgtacgtacgtacgtacgtacgtac")
})

test_that("parse_gb_records returns NA sequence when there is no ORIGIN block", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_true(is.na(parsed$sequence[parsed$accession == "GHI789"]))
})
