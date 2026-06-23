ncbi_results <- readRDS("results/ncbi_results.rds")

# Species with the most NCBI sequences (by accession count)
species_seq_count <- ncbi_results |>
  dplyr::filter(!is.na(organism), !is.na(accession)) |>
  dplyr::count(organism, name = "n_accessions") |>
  dplyr::arrange(dplyr::desc(n_accessions))

accessions_df <- ncbi_results |>
    dplyr::filter(organism == "Entomobrya nivalis") |>
    dplyr::pull(accession) |>
    unique()

# Fetch GenBank records (sequence + gene annotations) for each accession
library(rentrez)
set_entrez_key(Sys.getenv("NCBI_API_KEY"))

batch_size <- 50
batches <- split(accessions_df, ceiling(seq_along(accessions_df) / batch_size))

gb_records <- purrr::map(batches, \(batch) {
  rentrez::entrez_fetch(
    db       = "nucleotide",
    id       = batch,
    rettype  = "gb",
    retmode  = "text"
  )
}, .progress = TRUE)

# Parse gene name and sequence from each GenBank flat-file record
parse_gb_records <- function(gb_text) {
  # Split into individual records at the "//" delimiter
  records <- strsplit(gb_text, "(?m)^//\\s*$", perl = TRUE)[[1]]
  records <- records[nzchar(trimws(records))]

  purrr::map_dfr(records, \(rec) {
    accession <- regmatches(rec, regexpr("(?m)^ACCESSION\\s+(\\S+)", rec, perl = TRUE))
    accession <- if (length(accession)) sub("ACCESSION\\s+", "", accession) else NA_character_

    # Extract gene names from /gene="..." tags
    gene_matches <- regmatches(rec, gregexpr('/gene="([^"]+)"', rec, perl = TRUE))[[1]]
    genes <- unique(sub('/gene="([^"]+)"', "\\1", gene_matches))
    gene <- if (length(genes)) paste(genes, collapse = ";") else NA_character_

    # Extract the ORIGIN sequence block
    origin <- regmatches(rec, regexpr("(?ms)^ORIGIN.*", rec, perl = TRUE))
    sequence <- if (length(origin) && nzchar(origin)) {
      gsub("[^acgtACGT]", "", origin)
    } else {
      NA_character_
    }

    tibble::tibble(accession = accession, gene = gene, sequence = sequence)
  })
}

seq_data <- purrr::map_dfr(gb_records, parse_gb_records)

# Helper: check for in-frame stop codons in all 3 forward reading frames
has_stop_codon_coi <- function(seq) {
  if (is.na(seq)) return(NA)
  seq <- tolower(seq)
  codon_table <- Biostrings::getGeneticCode("SGC4")  # invertebrate mitochondrial code
  purrr::some(1:3, \(frame) {
    dna <- Biostrings::DNAString(seq)
    aa  <- Biostrings::translate(dna[frame:nchar(seq)], genetic.code = codon_table, if.fuzzy.codon = "solve")
    any(as.character(aa) == "*")
  })
}

# Sequence quality metrics
seq_data <- seq_data |>
  dplyr::mutate(
    seq_length     = nchar(sequence),
    n_count        = nchar(gsub("[^nN]", "", sequence)),
    n_pct          = round(n_count / seq_length * 100, 2),
    gc_pct         = round(
      (nchar(gsub("[^gcGC]", "", sequence)) / seq_length) * 100, 2
    )
  )

summary(seq_data$seq_length)

seq_data_filtered <- seq_data |>
  dplyr::filter(seq_length <= 680)

seqs <- DNAStringSet(seq_data_filtered$sequence)
aligns_seq <- DECIPHER::AlignSeqs(seqs)
BrowseSeqs(aligns_seq)
