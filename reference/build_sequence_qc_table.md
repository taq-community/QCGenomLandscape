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
  progress = TRUE,
  cores = NULL
)
```

## Arguments

- gb_records:

  List of raw GenBank flat-file blobs, e.g. from
  [`fetch_gb_records()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_gb_records.md)

- log_every_batch:

  Integer, chunk size for progress logging, parallel dispatch of batch
  parsing, AND per-chunk QC scoring (parse -\> score -\> stop-codon
  check all happen per chunk, not in three full-table passes), default
  500

- progress:

  Logical, also show `purrr` progress bars, default `TRUE`

- cores:

  Integer, cores to use for parallel batch parsing, default `NULL`
  meaning `min(4, parallel::detectCores() - 1)` (at least 1).
  Deliberately conservative rather than using every core: each
  `mclapply()` worker forks a copy of its batch text, so peak memory
  scales with the core count on top of the per-batch footprint –
  `detectCores() - 1` previously meant 15 concurrent forks on a 16-core
  machine, which is what pushed this step into the OOM killer under any
  other memory pressure on the same machine.

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
