#' Check for in-frame stop codons in a nucleotide sequence (COI / SGC4 code)
#'
#' Translates `seq` in all 3 forward reading frames using the invertebrate
#' mitochondrial genetic code (NCBI translation table 5 / SGC4, the standard
#' code for COI barcoding) and flags whether any frame contains a stop codon.
#'
#' @param seq Character scalar, a DNA sequence (case-insensitive)
#' @return Logical scalar, or `NA` if `seq` is `NA`
#' @examples
#' \dontrun{
#' has_stop_codon_coi("atgtaaatg")
#' }
#' @export
has_stop_codon_coi <- function(seq) {
  if (is.na(seq)) return(NA)
  seq <- tolower(seq)
  codon_table <- Biostrings::getGeneticCode("SGC4") # invertebrate mitochondrial code
  purrr::some(1:3, \(frame) {
    dna <- Biostrings::DNAString(seq)
    # a trailing incomplete codon (frame doesn't divide the sequence evenly)
    # is expected on almost every real sequence -- suppress that routine warning
    aa <- suppressWarnings(
      Biostrings::translate(dna[frame:nchar(seq)], genetic.code = codon_table, if.fuzzy.codon = "solve")
    )
    # as.character(aa) collapses the whole AAString into one string (e.g.
    # "K*K"), so a plain `== "*"` only matches a single-codon all-stop
    # sequence -- grepl() is needed to catch a stop codon anywhere in it
    grepl("*", as.character(aa), fixed = TRUE)
  })
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
