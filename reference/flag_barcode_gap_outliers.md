# Flag DNA sequences whose distance to conspecifics is a group outlier

The "barcode gap" idea from DNA barcoding QC: for each (species, group)
with at least `min_n` sequences, aligns them with
[`align_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/align_sequences.md)
(DECIPHER) and computes a pairwise distance matrix, then flags sequences
whose mean distance to the rest of the group is a robust (median + MAD)
outlier – a sequence far outside its own species' normal intraspecific
variation is a candidate misidentification or contamination, rather than
a fixed cross-taxon distance threshold, which wouldn't hold uniformly
across markers/taxa.

## Usage

``` r
flag_barcode_gap_outliers(
  sequence,
  species,
  group,
  k = 5,
  min_n = 3,
  max_n = 200,
  cores = NULL,
  log_every = 200
)
```

## Arguments

- sequence:

  Character vector of DNA sequences

- species, group:

  Vectors identifying species and marker/gene group; sequences are only
  compared within the same (species, group)

- k:

  Numeric, threshold in MAD units, default 5

- min_n:

  Integer, minimum (species, group) size to attempt, default 3

- max_n:

  Integer, (species, group) combinations larger than this are skipped
  rather than aligned, default 200

- cores:

  Integer, cores to use for parallel alignment, default `NULL` meaning
  `parallel::detectCores() - 1` (at least 1)

- log_every:

  Integer, log a progress line every this many aligned (species, group)
  combinations, default 200

## Value

Logical vector, same length as `sequence` – `NA` where not attempted
(group too small/large, alignment failed, or `NA` input)

## Details

This is the most expensive of the sequence-QC checks – multiple sequence
alignment scales worse than linearly with group size. Each (species,
group) aligns independently of every other, so the groups are
parallelized with
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
(fork-based – Linux/macOS only; falls back to a single core on Windows),
the same pattern used for batch parsing in
[`build_sequence_qc_table()`](https://taq-community.github.io/QCGenomLandscape/reference/build_sequence_qc_table.md);
`max_n` still skips groups too large to align cheaply rather than trying
and stalling one worker. Each worker's alignment is wrapped in
[`try()`](https://rdrr.io/r/base/try.html) – not just `mclapply()`'s own
error containment, which only isolates a failure when it actually forks
and silently stops containing errors at all if forking isn't available –
so one bad group can't take down the whole call.
