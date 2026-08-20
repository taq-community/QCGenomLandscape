# Flag sequence-length outliers within each marker/gene group

Robust (median + MAD) outlier detection, run independently per `group` –
a fixed length threshold isn't meaningful across markers as different as
a ~650bp COI barcode fragment and a ~1500bp rbcL sequence, and
hardcoding a per-marker "expected length" table would need constant
upkeep as new marker types show up in the data. MAD-based (not SD-based)
specifically because SD is itself inflated by the outliers being
detected; groups smaller than `min_n` are left unflagged (`NA`) since
there isn't enough data in a tiny group to say what's "normal".

## Usage

``` r
flag_length_outliers(seq_length, group, k = 5, min_n = 10)
```

## Arguments

- seq_length:

  Numeric vector of sequence lengths

- group:

  Vector (e.g. `gene_group` from
  [`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md))
  identifying which group each length belongs to

- k:

  Numeric, threshold in MAD units, default 5 (a conservative bar – under
  normality this flags roughly the most extreme ~0.0001% of a
  distribution, well past ordinary biological length variation)

- min_n:

  Integer, minimum group size to attempt outlier detection, default 10

## Value

Logical vector, same length as `seq_length` – `TRUE` for outliers, `NA`
where `seq_length` is `NA` or the group is smaller than `min_n`

## Examples

``` r
flag_length_outliers(c(640, 655, 648, 651, 5), rep("COI", 5), min_n = 3)
#> [1] FALSE FALSE FALSE FALSE  TRUE
```
