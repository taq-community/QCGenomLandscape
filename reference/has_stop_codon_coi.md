# Flag DNA sequences with no clean reading frame (COI / SGC4 code)

Translates each sequence in all 3 forward reading frames using the
invertebrate mitochondrial genetic code (NCBI translation table 5 /
SGC4, the standard code for COI barcoding), and flags a sequence when
**none** of the 3 frames translates without hitting a stop codon.

## Usage

``` r
has_stop_codon_coi(seq)
```

## Arguments

- seq:

  Character vector of DNA sequences (case-insensitive); `NA` elements
  return `NA`

## Value

Logical vector, the same length as `seq` – `TRUE` means no reading frame
translates cleanly (likely problem)

## Details

This is deliberately *not* "does any frame contain a stop" – for a real
~650bp protein-coding sequence, a random frame is stop-free with
probability roughly `(61/64)^~200`, i.e. astronomically small by chance,
so at least one clean frame exists for essentially every real,
correctly-oriented sequence. Flagging "a stop exists in some frame" is
true of ~100% of real sequences checked this way (3 chances for an
unrelated frame to hit one of 4 stop codons) and isn't discriminating.
Flagging "no frame is clean" is the rare, meaningful signal instead –
consistent with a frameshift, NUMT/pseudogene, wrong orientation, or
sequencing error.

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
#> [1] FALSE
has_stop_codon_coi(c("atgtaaatg", "aaaaaaaaa", NA))
#> [1] FALSE FALSE    NA
```
