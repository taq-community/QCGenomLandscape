# Summarise the species lost to failed Entrez queries

[`fetch_ncbi_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_ncbi_sequences.md)
records every query it could not complete in `deficient_queries` and
carries on, so a failed batch does not abort a multi-hour run. Nothing
downstream consumed that record, which is how a single timed-out batch
became a published "Amphibians: 0% coverage" panel. This turns the
failure log back into the list of species it cost, so the pipeline can
surface – or refuse to pass – an incomplete corpus.

## Usage

``` r
deficient_query_species(deficient_queries, species)
```

## Arguments

- deficient_queries:

  List as returned in `fetch_ncbi_sequences()$deficient_queries`

- species:

  Character vector of species names to look for in the failed query
  strings, typically the BDQC species list

## Value

Tibble with one row per (query_index, species) lost, columns
`query_index`, `error_type`, `error_message`, `species`
