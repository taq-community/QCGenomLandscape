# Align a set of DNA sequences

Thin, namespaced wrapper over
[`Biostrings::DNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html) +
[`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html).

## Usage

``` r
align_sequences(sequences, verbose = FALSE)
```

## Arguments

- sequences:

  Character vector of DNA sequences

- verbose:

  Logical, print
  [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html)'s
  progress output, default `FALSE` (it's chatty by default – fine for
  interactive use, unwanted noise when called many times in a loop, e.g.
  from
  [`flag_barcode_gap_outliers()`](https://taq-community.github.io/QCGenomLandscape/reference/flag_barcode_gap_outliers.md))

## Value

A
[`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
of aligned sequences
