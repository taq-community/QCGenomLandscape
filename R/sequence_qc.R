#' Check for in-frame stop codons in DNA sequences (COI / SGC4 code)
#'
#' Translates each sequence in all 3 forward reading frames using the
#' invertebrate mitochondrial genetic code (NCBI translation table 5 / SGC4,
#' the standard code for COI barcoding) and flags whether any frame contains
#' a stop codon.
#'
#' Vectorized over `seq`: all sequences are translated together per frame in
#' one `Biostrings::translate()` call (3 calls total) instead of looping
#' element by element and re-deriving the genetic code table each time. This
#' is the hot path in [build_sequence_qc_table()] -- a full BDQC run checks
#' hundreds of thousands of sequences, where the old per-element scalar
#' version took hours.
#'
#' @param seq Character vector of DNA sequences (case-insensitive); `NA`
#'   elements return `NA`
#' @return Logical vector, the same length as `seq`
#' @examples
#' has_stop_codon_coi("atgtaaatg")
#' has_stop_codon_coi(c("atgtaaatg", "aaaaaaaaa", NA))
#' @export
has_stop_codon_coi <- function(seq) {
  result <- rep(NA, length(seq))
  valid <- !is.na(seq)
  if (!any(valid)) {
    return(result)
  }

  dna_set <- Biostrings::DNAStringSet(tolower(seq[valid]))
  codon_table <- Biostrings::getGeneticCode("SGC4") # invertebrate mitochondrial code
  seq_width <- Biostrings::width(dna_set)

  has_stop <- logical(length(dna_set))
  for (frame in 1:3) {
    # sequences shorter than the frame offset have no codons in that frame
    in_frame <- seq_width >= frame
    if (!any(in_frame)) next

    sub_set <- Biostrings::subseq(dna_set[in_frame], start = frame)
    # a trailing incomplete codon (frame doesn't divide the sequence evenly)
    # is expected on almost every real sequence -- suppress that routine warning
    aa <- suppressWarnings(
      Biostrings::translate(sub_set, genetic.code = codon_table, if.fuzzy.codon = "solve")
    )
    # as.character(aa) collapses each AAString into one string (e.g. "K*K"),
    # so a plain `== "*"` only matches a single-codon all-stop sequence --
    # grepl() is needed to catch a stop codon anywhere in it
    has_stop[in_frame] <- has_stop[in_frame] | grepl("*", as.character(aa), fixed = TRUE)
  }

  result[valid] <- has_stop
  result
}

#' Compute basic sequence-quality metrics
#'
#' @param sequence Character vector of DNA sequences
#' @return Tibble with columns `seq_length`, `n_count`, `n_pct`, `gc_pct`
#'   (one row per input sequence)
#' @examples
#' score_sequence_quality(c("acgtacgt", "acgtnnnn"))
#' @export
score_sequence_quality <- function(sequence) {
  seq_length <- nchar(sequence)
  n_count <- nchar(gsub("[^nN]", "", sequence))
  n_pct <- round(n_count / seq_length * 100, 2)
  gc_pct <- round(nchar(gsub("[^gcGC]", "", sequence)) / seq_length * 100, 2)

  tibble::tibble(
    seq_length = seq_length,
    n_count = n_count,
    n_pct = n_pct,
    gc_pct = gc_pct
  )
}

#' Build the per-record sequence-QC table from fetched GenBank batches
#'
#' Parses raw GenBank flat-file batches, scores basic sequence-quality
#' metrics, and flags in-frame stop codons (COI/SGC4 code). Logs progress
#' periodically via `logger` -- at BDQC scale (thousands of batches,
#' hundreds of thousands of sequences) this step has no other progress
#' signal and can run for hours, particularly the per-sequence stop-codon
#' check (one `Biostrings` translation call per sequence).
#'
#' Batch parsing is parallelized with [parallel::mclapply()] (fork-based --
#' Linux/macOS only; falls back to a single core on Windows, since forking
#' isn't available there) since each batch parses independently of the
#' others.
#'
#' @param gb_records List of raw GenBank flat-file blobs, e.g. from
#'   [fetch_gb_records()]
#' @param log_every_batch Integer, chunk size for both progress logging and
#'   parallel dispatch of the batch-parsing step, default 500
#' @param log_every_seq Integer, log a progress line every this many
#'   stop-codon checks, default 20000
#' @param progress Logical, also show `purrr` progress bars, default `TRUE`
#' @param cores Integer, cores to use for parallel batch parsing, default
#'   `NULL` meaning `parallel::detectCores() - 1` (at least 1). Deliberately
#'   conservative rather than using every core, since each worker holds its
#'   own copy of the batch text and this step has previously run under heavy
#'   memory pressure.
#' @return Tibble with columns `accession`, `gene`, `sequence`, `seq_length`,
#'   `n_count`, `n_pct`, `gc_pct`, `has_stop`
#' @export
build_sequence_qc_table <- function(gb_records,
                                     log_every_batch = 500,
                                     log_every_seq = 20000,
                                     progress = TRUE,
                                     cores = NULL) {
  n_batches <- length(gb_records)
  if (is.null(cores)) {
    detected <- parallel::detectCores()
    cores <- if (is.na(detected)) 1L else max(1L, detected - 1L)
  }
  if (.Platform$OS.type == "windows" && cores > 1) {
    logger::log_info("seq_data: parallel batch parsing needs fork(), unavailable on Windows -- using 1 core")
    cores <- 1
  }
  logger::log_info("seq_data: parsing {n_batches} GenBank batches ({cores} core(s))...")

  batch_chunks <- split(seq_len(n_batches), ceiling(seq_len(n_batches) / log_every_batch))
  parsed <- purrr::map_dfr(batch_chunks, function(idx) {
    # try() (not just mclapply()'s own error containment) since mclapply()
    # only isolates errors when it actually forks -- it silently falls back
    # to plain lapply() (no containment at all) when forking isn't
    # available, so relying on it alone would let one bad batch kill the run
    chunk_results <- parallel::mclapply(gb_records[idx], function(batch) {
      try(parse_gb_records(batch), silent = TRUE)
    }, mc.cores = cores)

    failed <- vapply(chunk_results, function(x) inherits(x, "try-error"), logical(1))
    if (any(failed)) {
      logger::log_warn(
        "seq_data: {sum(failed)} batch(es) up to {max(idx)}/{n_batches} failed to parse -- dropping"
      )
      chunk_results <- chunk_results[!failed]
    }

    logger::log_info("seq_data: parsed batch {max(idx)}/{n_batches}")
    dplyr::bind_rows(chunk_results)
  }, .progress = progress)

  logger::log_info("seq_data: parsed {nrow(parsed)} records -- scoring sequence quality...")
  qc <- dplyr::bind_cols(parsed, score_sequence_quality(parsed$sequence))

  n_seq <- nrow(qc)
  logger::log_info("seq_data: checking {n_seq} sequences for in-frame stop codons...")
  # has_stop_codon_coi() is itself vectorized (one Biostrings::translate()
  # call per frame, not per sequence) -- chunking here is only to keep
  # progress visible on a full run, not for speed
  chunks <- split(seq_len(n_seq), ceiling(seq_len(n_seq) / log_every_seq))
  has_stop <- vector("logical", n_seq)
  for (idx in chunks) {
    has_stop[idx] <- has_stop_codon_coi(qc$sequence[idx])
    logger::log_info("seq_data: stop-codon check {max(idx)}/{n_seq}")
  }

  logger::log_success("seq_data: done ({n_seq} records)")
  dplyr::mutate(qc, has_stop = has_stop)
}

#' Align a set of DNA sequences
#'
#' Thin, namespaced wrapper over `Biostrings::DNAStringSet()` +
#' `DECIPHER::AlignSeqs()`.
#'
#' @param sequences Character vector of DNA sequences
#' @return A `Biostrings::DNAStringSet` of aligned sequences
#' @export
align_sequences <- function(sequences) {
  DECIPHER::AlignSeqs(Biostrings::DNAStringSet(sequences))
}
