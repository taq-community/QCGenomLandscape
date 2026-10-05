# Extract a `Genus species` binomial, preserving vector length

[`regmatches()`](https://rdrr.io/r/base/regmatches.html) drops
non-matching elements, so using it inside
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
turns one malformed name into a whole-pipeline size error. This returns
`NA` for rows that carry no binomial instead.

## Usage

``` r
extract_binomial(x)
```

## Arguments

- x:

  Character vector of scientific names

## Value

Character vector the same length as `x`

## Examples

``` r
extract_binomial(c("Ambystoma maculatum", "Ambystoma (2) laterale - jeffersonianum"))
#> [1] "Ambystoma maculatum" NA                   
```
