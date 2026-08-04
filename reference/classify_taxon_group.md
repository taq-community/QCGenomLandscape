# Classify BDQC taxonomic groups into finer display categories

Reclassifies the coarse `group_en` column of the BDQC species list into
finer taxonomic categories used for figures (e.g. splitting mammals into
marine/terrestrial, arthropods into insects/crustaceans).

## Usage

``` r
classify_taxon_group(
  group_en,
  order = NA,
  class = NA,
  phylum = NA,
  kingdom = NA
)
```

## Arguments

- group_en:

  Character vector, BDQC's coarse taxonomic group

- order:

  Character vector, taxonomic order (same length as `group_en`)

- class:

  Character vector, taxonomic class (same length as `group_en`)

- phylum:

  Character vector, taxonomic phylum (same length as `group_en`)

- kingdom:

  Character vector, taxonomic kingdom (same length as `group_en`)

## Value

Character vector of finer taxonomic group labels

## Examples

``` r
classify_taxon_group("Mammals", order = "Cetacea")
#> [1] "Marine mammals"
classify_taxon_group("Arthropods", class = "Insecta")
#> [1] "Insects"
```
