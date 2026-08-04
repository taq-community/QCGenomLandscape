# Compute basic sequence-quality metrics

Compute basic sequence-quality metrics

## Usage

``` r
score_sequence_quality(sequence)
```

## Arguments

- sequence:

  Character vector of DNA sequences

## Value

Tibble with columns `seq_length`, `n_count`, `n_pct`, `gc_pct` (one row
per input sequence)

## Examples

``` r
score_sequence_quality(c("acgtacgt", "acgtnnnn"))
#> # A tibble: 2 × 4
#>   seq_length n_count n_pct gc_pct
#>        <int>   <int> <dbl>  <dbl>
#> 1          8       0     0     50
#> 2          8       4    50     25
```
