# Donut chart: genomic-data coverage by species-at-risk status

Donut chart: genomic-data coverage by species-at-risk status

## Usage

``` r
plot_risk_status_coverage(risk_genes, jurisdiction = "QC")
```

## Arguments

- risk_genes:

  Tibble with columns `species`, `status`, `jurisdiction`, `has_data`
  (risk status joined to gene-data presence, e.g. via
  [`load_risk_status()`](https://taq-community.github.io/QCGenomLandscape/reference/load_risk_status.md)
  joined with a `gene_group` presence flag)

- jurisdiction:

  Character, which `jurisdiction` value to filter to, default `"QC"`

## Value

A `ggplot` object (caller is responsible for
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html))
