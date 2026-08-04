# Align a set of DNA sequences

Thin, namespaced wrapper over
[`Biostrings::DNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html) +
[`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html).

## Usage

``` r
align_sequences(sequences)
```

## Arguments

- sequences:

  Character vector of DNA sequences

## Value

A
[`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
of aligned sequences
