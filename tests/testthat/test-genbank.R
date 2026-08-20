test_that("parse_gb_records parses accession, gene, and sequence per record", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_equal(nrow(parsed), 3)
  expect_equal(parsed$accession, c("ABC123", "DEF456", "GHI789"))
})

test_that("parse_gb_records extracts a single-line DEFINITION", {
  gb_text <- paste(
    readLines(testthat::test_path("fixtures", "sample_genbank_record.gb")),
    collapse = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_equal(parsed$definition, c(
    "Test species mitochondrial COI gene, single gene tag.",
    "Test species with two distinct gene tags.",
    "Test species without a sequence (no ORIGIN block)."
  ))
})

test_that("parse_gb_records unwraps a DEFINITION that spans indented continuation lines", {
  gb_text <- paste(
    "LOCUS       XYZ999               16746 bp    DNA     circular VRT 18-DEC-2019",
    "DEFINITION  Boreogadus saida voucher UW151986-5 mitochondrion,",
    "            complete genome.",
    "ACCESSION   XYZ999",
    "VERSION     XYZ999.1",
    "FEATURES             Location/Qualifiers",
    "     gene            1..10",
    "                     /gene=\"COI\"",
    "ORIGIN",
    "        1 acgtacgtac",
    "//",
    sep = "\n"
  )

  parsed <- parse_gb_records(gb_text)

  expect_equal(parsed$definition, "Boreogadus saida voucher UW151986-5 mitochondrion, complete genome.")
})

test_that("is_complete_genome matches GenBank's own 'complete genome' wording, not 'partial genome'", {
  out <- is_complete_genome(c(
    "Boreogadus saida mitochondrion, complete genome.",
    "Boreogadus saida mitochondrion, partial genome.",
    "Boreogadus saida cytochrome oxidase subunit I (COI) gene, partial cds.",
    NA
  ))

  expect_equal(out, c(TRUE, FALSE, FALSE, NA))
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

# ---- fetch_gb_records: stubbed fetch_fn, no network ----

test_that("fetch_gb_records returns one blob per batch when every fetch succeeds", {
  out <- fetch_gb_records(
    accessions = c("A1", "A2", "A3"),
    batch_size = 2,
    fetch_fn = function(db, id, rettype, retmode) paste(id, collapse = ","),
    progress = FALSE
  )

  expect_length(out, 2)
  expect_equal(out[[1]], "A1,A2")
  expect_equal(out[[2]], "A3")
})

test_that("fetch_gb_records retries a failing batch and keeps the result once it succeeds", {
  attempts <- 0
  flaky_fetch_fn <- function(db, id, rettype, retmode) {
    attempts <<- attempts + 1
    if (attempts < 2) stop("simulated transient SSL error")
    "recovered text"
  }

  out <- fetch_gb_records(
    accessions = c("A1"),
    batch_size = 50,
    fetch_fn = flaky_fetch_fn,
    max_retries = 3,
    sleep_fn = function(seconds) invisible(NULL),
    progress = FALSE
  )

  expect_equal(attempts, 2)
  expect_equal(out, list("recovered text"))
})

test_that("fetch_gb_records drops a batch that fails on every attempt instead of erroring", {
  always_fails_fn <- function(db, id, rettype, retmode) stop("simulated persistent SSL error")

  out <- fetch_gb_records(
    accessions = c("A1", "A2", "A3", "A4"),
    batch_size = 2,
    fetch_fn = function(db, id, rettype, retmode) {
      if (identical(id, c("A1", "A2"))) always_fails_fn(db, id, rettype, retmode)
      paste(id, collapse = ",")
    },
    max_retries = 2,
    sleep_fn = function(seconds) invisible(NULL),
    progress = FALSE
  )

  expect_length(out, 1)
  expect_equal(out[[1]], "A3,A4")
})
