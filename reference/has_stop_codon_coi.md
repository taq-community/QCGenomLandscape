# Check for in-frame stop codons in DNA sequences (COI / SGC4 code)

Translates each sequence in all 3 forward reading frames using the
invertebrate mitochondrial genetic code (NCBI translation table 5 /
SGC4, the standard code for COI barcoding) and flags whether any frame
contains a stop codon.

## Usage

``` r
has_stop_codon_coi(seq)
```

## Arguments

- seq:

  Character vector of DNA sequences (case-insensitive); `NA` elements
  return `NA`

## Value

Logical vector, the same length as `seq`

## Details

Vectorized over `seq`: all sequences are translated together per frame
in one
[`Biostrings::translate()`](https://rdrr.io/pkg/Biostrings/man/translate.html)
call (3 calls total) instead of looping element by element and
re-deriving the genetic code table each time. This is the hot path in
[`build_sequence_qc_table()`](https://taq-community.github.io/QCGenomLandscape/reference/build_sequence_qc_table.md)
– a full BDQC run checks hundreds of thousands of sequences, where the
old per-element scalar version took hours.

## Examples

``` r
has_stop_codon_coi("atgtaaatg")
#> [1] TRUE
has_stop_codon_coi(c("atgtaaatg", "aaaaaaaaa", NA))
#> [1]  TRUE FALSE    NA
```
