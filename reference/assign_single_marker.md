# Resolve a gene annotation to a single comparable marker

[`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md)
is deliberately coarse: three of its labels are *buckets* of several
distinct markers – `"Photosynthesis-related"` pools rbcL (~1400 bp),
matK (~850 bp) and trnL (~300 bp), `"Nuclear rRNA / ITS"` pools 5.8S
(~160 bp) with LSU (~3.5 kb), and `"Fungal protein-coding"` pools RPB1
with RPB2. That is fine for coverage figures ("does this species have a
plastid marker at all?") but wrong for any length-based quality metric,
which needs sequences that are actually comparable to each other.

## Usage

``` r
assign_single_marker(gene, definition = NULL)
```

## Arguments

- gene:

  Character vector of raw gene names, as annotated on the record

- definition:

  Character vector of GenBank DEFINITION lines, same length as `gene`;
  passed through to
  [`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md)
  so that records with no `/gene=` tag (overwhelmingly ITS) resolve to
  `"ITS"` instead of dropping out of every per-marker statistic. Default
  `NULL`.

## Value

Character vector of marker labels, `NA` where no single marker applies

## Details

This resolves those buckets down to the raw gene name (normalised for
case and for the optional ` rRNA` suffix, so `"18S rRNA"` and `"18S"`
are one marker), and keeps
[`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md)'s
label everywhere else – that is where the label already denotes exactly
one marker *and* usefully unifies synonyms (`COX1`/`COI`, `cytb`/`cob`).

Genome-scale records and the `"Other"` catch-all resolve to `NA`:
neither is a single marker, so neither belongs in a per-marker
comparison.

## Examples

``` r
assign_single_marker(c("COX1", "rbcL", "matK", "18S rRNA", "18S"))
#> [1] "COI"  "RBCL" "MATK" "18S"  "18S" 
assign_single_marker("ND1;ND2;COX1") # NA -- genome-scale record
#> [1] NA
```
