#' Flag DNA sequences with no clean reading frame (COI / SGC4 code)
#'
#' Translates each sequence in all 3 forward reading frames using the
#' invertebrate mitochondrial genetic code (NCBI translation table 5 / SGC4,
#' the standard code for COI barcoding), and flags a sequence when **none**
#' of the 3 frames translates without hitting a stop codon.
#'
#' This is deliberately *not* "does any frame contain a stop" -- for a real
#' ~650bp protein-coding sequence, a random frame is stop-free with
#' probability roughly `(61/64)^~200`, i.e. astronomically small by chance,
#' so at least one clean frame exists for essentially every real,
#' correctly-oriented sequence. Flagging "a stop exists in some frame" is
#' true of ~100% of real sequences checked this way (3 chances for an
#' unrelated frame to hit one of 4 stop codons) and isn't discriminating.
#' Flagging "no frame is clean" is the rare, meaningful signal instead --
#' consistent with a frameshift, NUMT/pseudogene, wrong orientation, or
#' sequencing error.
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
#' @return Logical vector, the same length as `seq` -- `TRUE` means no
#'   reading frame translates cleanly (likely problem)
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

  # frame_clean[, f] == TRUE means frame f translates with no stop codon
  frame_clean <- matrix(FALSE, nrow = length(dna_set), ncol = 3)
  for (frame in 1:3) {
    # sequences shorter than the frame offset have no codons in that frame
    # -- leave them FALSE (not clean) rather than claim a frame that doesn't exist
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
    frame_clean[in_frame, frame] <- !grepl("*", as.character(aa), fixed = TRUE)
  }

  has_clean_frame <- frame_clean[, 1] | frame_clean[, 2] | frame_clean[, 3]
  result[valid] <- !has_clean_frame
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

#' Flag sequence-length outliers within each marker/gene group
#'
#' Robust (median + MAD) outlier detection, run independently per `group` --
#' a fixed length threshold isn't meaningful across markers as different as
#' a ~650bp COI barcode fragment and a ~1500bp rbcL sequence, and hardcoding
#' a per-marker "expected length" table would need constant upkeep as new
#' marker types show up in the data. MAD-based (not SD-based) specifically
#' because SD is itself inflated by the outliers being detected; groups
#' smaller than `min_n` are left unflagged (`NA`) since there isn't enough
#' data in a tiny group to say what's "normal".
#'
#' @param seq_length Numeric vector of sequence lengths
#' @param group Vector (e.g. `gene_group` from [assign_gene_group()])
#'   identifying which group each length belongs to
#' @param k Numeric, threshold in MAD units, default 5 (a conservative bar --
#'   under normality this flags roughly the most extreme ~0.0001% of a
#'   distribution, well past ordinary biological length variation)
#' @param min_n Integer, minimum group size to attempt outlier detection,
#'   default 10
#' @return Logical vector, same length as `seq_length` -- `TRUE` for
#'   outliers, `NA` where `seq_length` is `NA` or the group is smaller than
#'   `min_n`
#' @examples
#' flag_length_outliers(c(640, 655, 648, 651, 5), rep("COI", 5), min_n = 3)
#' @export
flag_length_outliers <- function(seq_length, group, k = 5, min_n = 10) {
  stopifnot(length(seq_length) == length(group))

  result <- rep(NA, length(seq_length))
  idx_by_group <- split(seq_along(seq_length), group)

  for (idx in idx_by_group) {
    valid <- idx[!is.na(seq_length[idx])]
    if (length(valid) < min_n) next

    lengths <- seq_length[valid]
    med <- stats::median(lengths)
    mad_val <- stats::mad(lengths)

    result[valid] <- if (mad_val == 0) {
      # degenerate: every length in the group is (near-)identical --
      # fall back to exact-match instead of dividing by zero
      lengths != med
    } else {
      abs(lengths - med) / mad_val > k
    }
  }

  result
}

#' Flag DNA sequences whose distance to conspecifics is a group outlier
#'
#' The "barcode gap" idea from DNA barcoding QC: for each (species, group)
#' with at least `min_n` sequences, aligns them with [align_sequences()]
#' (DECIPHER) and computes a pairwise distance matrix, then flags sequences
#' whose mean distance to the rest of the group is a robust (median + MAD)
#' outlier -- a sequence far outside its own species' normal intraspecific
#' variation is a candidate misidentification or contamination, rather than
#' a fixed cross-taxon distance threshold, which wouldn't hold uniformly
#' across markers/taxa.
#'
#' This is the most expensive of the sequence-QC checks -- multiple sequence
#' alignment scales worse than linearly with group size. Each (species,
#' group) aligns independently of every other, so the groups are
#' parallelized with [parallel::mclapply()] (fork-based -- Linux/macOS
#' only; falls back to a single core on Windows), the same pattern used for
#' batch parsing in [build_sequence_qc_table()]; `max_n` still skips groups
#' too large to align cheaply rather than trying and stalling one worker.
#' Each worker's alignment is wrapped in `try()` -- not just `mclapply()`'s
#' own error containment, which only isolates a failure when it actually
#' forks and silently stops containing errors at all if forking isn't
#' available -- so one bad group can't take down the whole call.
#'
#' @param sequence Character vector of DNA sequences
#' @param species,group Vectors identifying species and marker/gene group;
#'   sequences are only compared within the same (species, group)
#' @param k Numeric, threshold in MAD units, default 5
#' @param min_n Integer, minimum (species, group) size to attempt, default 3
#' @param max_n Integer, (species, group) combinations larger than this are
#'   skipped rather than aligned, default 200
#' @param cores Integer, cores to use for parallel alignment, default
#'   `NULL` meaning `parallel::detectCores() - 1` (at least 1)
#' @param log_every Integer, log a progress line every this many aligned
#'   (species, group) combinations, default 200
#' @return Logical vector, same length as `sequence` -- `NA` where not
#'   attempted (group too small/large, alignment failed, or `NA` input)
#' @export
flag_barcode_gap_outliers <- function(sequence, species, group, k = 5, min_n = 3, max_n = 200,
                                       cores = NULL, log_every = 200) {
  stopifnot(length(sequence) == length(species), length(sequence) == length(group))

  if (is.null(cores)) {
    detected <- parallel::detectCores()
    cores <- if (is.na(detected)) 1L else max(1L, detected - 1L)
  }
  if (.Platform$OS.type == "windows" && cores > 1) {
    logger::log_info("flag_barcode_gap_outliers: parallel alignment needs fork(), unavailable on Windows -- using 1 core")
    cores <- 1
  }

  key <- paste(species, group, sep = "\u0001")
  idx_by_key <- split(seq_along(sequence), key)
  # only attempt groups within [min_n, max_n] -- skip the rest up front so
  # they don't cost a worker dispatch just to be discarded
  eligible <- Filter(function(idx) {
    n_valid <- sum(!is.na(sequence[idx]))
    n_valid >= min_n && n_valid <= max_n
  }, idx_by_key)

  n_eligible <- length(eligible)
  logger::log_info(
    "flag_barcode_gap_outliers: aligning {n_eligible}/{length(idx_by_key)} ",
    "(species, group) combinations ({cores} core(s))..."
  )

  align_one <- function(idx) {
    try(
      {
        valid <- idx[!is.na(sequence[idx])]
        aligned <- align_sequences(sequence[valid])
        dist_mat <- DECIPHER::DistanceMatrix(aligned, verbose = FALSE, processors = 1)
        mean_dist <- rowMeans(dist_mat, na.rm = TRUE)
        med <- stats::median(mean_dist)
        mad_val <- stats::mad(mean_dist)
        flags <- if (mad_val == 0) rep(FALSE, length(mean_dist)) else (abs(mean_dist - med) / mad_val) > k
        stats::setNames(flags, valid)
      },
      silent = TRUE
    )
  }

  # chunked (not one big mclapply()) purely so progress is visible on a
  # long run -- mclapply() itself has no progress-reporting hook
  chunks <- split(eligible, ceiling(seq_along(eligible) / log_every))
  flags_by_group <- list()
  n_done <- 0
  for (chunk in chunks) {
    flags_by_group <- c(flags_by_group, parallel::mclapply(chunk, align_one, mc.cores = cores))
    n_done <- n_done + length(chunk)
    logger::log_info("flag_barcode_gap_outliers: aligned {n_done}/{n_eligible}")
  }

  result <- rep(NA, length(sequence))
  n_failed <- 0
  for (flags in flags_by_group) {
    if (inherits(flags, "try-error")) {
      n_failed <- n_failed + 1
      next
    }
    result[as.integer(names(flags))] <- unname(flags)
  }
  if (n_failed > 0) {
    logger::log_warn("flag_barcode_gap_outliers: {n_failed} group(s) failed to align -- left as NA")
  }

  result
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
#' @return Tibble with columns `accession`, `definition`, `gene`, `sequence`,
#'   `seq_length`, `n_count`, `n_pct`, `gc_pct`, `has_stop`,
#'   `is_complete_genome` (see [is_complete_genome()] -- GenBank's own
#'   "complete genome" wording, not inferred from the gene count)
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
  qc <- dplyr::bind_cols(parsed, score_sequence_quality(parsed$sequence)) |>
    dplyr::mutate(is_complete_genome = is_complete_genome(definition))

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
#' @param verbose Logical, print `DECIPHER::AlignSeqs()`'s progress output,
#'   default `FALSE` (it's chatty by default -- fine for interactive use,
#'   unwanted noise when called many times in a loop, e.g. from
#'   [flag_barcode_gap_outliers()])
#' @return A `Biostrings::DNAStringSet` of aligned sequences
#' @export
align_sequences <- function(sequences, verbose = FALSE) {
  DECIPHER::AlignSeqs(Biostrings::DNAStringSet(sequences), verbose = verbose)
}
