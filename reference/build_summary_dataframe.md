# Build the per-species genomic-data summary table

Pure, in-memory core of `create_dataframe.R`: joins NCBI sequence
counts, gene-marker counts, BDQC taxonomy, and conservation-risk status
into one row per species.

## Usage

``` r
build_summary_dataframe(ncbi_results, genes_df, bdqc_taxo, ca_risk, qc_risk)
```

## Arguments

- ncbi_results:

  Tibble as returned by `fetch_ncbi_sequences()$results` (must have
  `organism`, `accession` columns)

- genes_df:

  Tibble as returned by
  [`fetch_gene_annotations()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_gene_annotations.md)
  (must have `accession`, `gene` columns)

- bdqc_taxo:

  Tibble with columns `species`, `vernacular_fr`, `vernacular_en`,
  `group_en`

- ca_risk, qc_risk:

  Tibble with columns `species`, `status` (e.g. from
  [`load_risk_status()`](https://taq-community.github.io/QCGenomLandscape/reference/load_risk_status.md))

## Value

Tibble, one row per species
