# Donut chart: gene prevalence by taxonomic group

Donut chart: gene prevalence by taxonomic group

## Usage

``` r
plot_gene_prevalence(
  genes_grouped,
  gene_groups = c("COI", "Photosynthesis-related (rbcL, matK, etc.)"),
  min_species = 5
)
```

## Arguments

- genes_grouped:

  Tibble with columns `species`, `gene_group`, `groupe` (gene
  classification joined to species and their finer taxonomic group, e.g.
  via
  [`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md)
  and
  [`classify_taxon_group()`](https://taq-community.github.io/QCGenomLandscape/reference/classify_taxon_group.md))

- gene_groups:

  Character vector of `gene_group` values to plot as separate facets,
  default `c("COI", "Photosynthesis-related (rbcL, matK, etc.)")`

- min_species:

  Integer; gene groups with fewer than this many total species (summed
  across taxonomic groups) are dropped, default 5

## Value

A `ggplot` object (caller is responsible for
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html))
