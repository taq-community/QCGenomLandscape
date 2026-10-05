# Classify a GenBank DEFINITION line into a coarse marker group

Fallback for records with no `/gene=` qualifier, used by
[`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md).
Order matters: a single ITS deposit routinely names the small subunit,
ITS1, 5.8S, ITS2 and the large subunit in one DEFINITION, so the
nuclear-rRNA test has to come before the narrower ribosomal tests or
those records scatter across several groups.

## Usage

``` r
classify_definition_group(definition)
```

## Arguments

- definition:

  Character vector of GenBank DEFINITION lines

## Value

Character vector of
[`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md)
labels

## Examples

``` r
classify_definition_group(
  "Laccaria amethystina internal transcribed spacer 1, partial sequence"
)
#> [1] "Nuclear rRNA / ITS"
```
