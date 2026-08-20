# Build the per-record sequence-QC table from fetched GenBank batches

Parses raw GenBank flat-file batches, scores basic sequence-quality
metrics, and flags in-frame stop codons (COI/SGC4 code). Logs progress
periodically via `logger` – at BDQC scale (thousands of batches,
hundreds of thousands of sequences) this step has no other progress
signal and can run for hours, particularly the per-sequence stop-codon
check (one `Biostrings` translation call per sequence).

## Usage

``` r
build_sequence_qc_table(
  gb_records,
  log_every_batch = 500,
  log_every_seq = 20000,
  progress = TRUE,
  cores = NULL
)
```

## Arguments

- gb_records:

  List of raw GenBank flat-file blobs, e.g. from
  [`fetch_gb_records()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_gb_records.md)

- log_every_batch:

  Integer, chunk size for both progress logging and parallel dispatch of
  the batch-parsing step, default 500

- log_every_seq:

  Integer, log a progress line every this many stop-codon checks,
  default 20000

- progress:

  Logical, also show `purrr` progress bars, default `TRUE`

- cores:

  Integer, cores to use for parallel batch parsing, default `NULL`
  meaning `parallel::detectCores() - 1` (at least 1). Deliberately
  conservative rather than using every core, since each worker holds its
  own copy of the batch text and this step has previously run under
  heavy memory pressure.

## Value

Tibble with columns `accession`, `definition`, `gene`, `sequence`,
`seq_length`, `n_count`, `n_pct`, `gc_pct`, `has_stop`,
`is_complete_genome` (see
[`is_complete_genome()`](https://taq-community.github.io/QCGenomLandscape/reference/is_complete_genome.md)
– GenBank's own "complete genome" wording, not inferred from the gene
count)

## Details

Batch parsing is parallelized with
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
(fork-based – Linux/macOS only; falls back to a single core on Windows,
since forking isn't available there) since each batch parses
independently of the others.
