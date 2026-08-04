# Check for in-frame stop codons in a nucleotide sequence (COI / SGC4 code)

Translates `seq` in all 3 forward reading frames using the invertebrate
mitochondrial genetic code (NCBI translation table 5 / SGC4, the
standard code for COI barcoding) and flags whether any frame contains a
stop codon.

## Usage

``` r
has_stop_codon_coi(seq)
```

## Arguments

- seq:

  Character scalar, a DNA sequence (case-insensitive)

## Value

Logical scalar, or `NA` if `seq` is `NA`

## Examples

``` r
if (FALSE) { # \dontrun{
has_stop_codon_coi("atgtaaatg")
} # }
```
