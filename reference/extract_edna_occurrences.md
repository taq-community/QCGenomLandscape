# Normalize NCBI + BOLD records into a common eDNA occurrence table

Pulls the fields needed to compare eDNA records against traditional
occurrence data (e.g. the Biodiversité Québec Atlas) spatio-temporally:
one row per geo-referenced eDNA record, with species, coordinates, and
collection date normalized to a common shape regardless of source.
Records with no parseable coordinates, or no species name, are dropped –
both sources have some (BOLD: unidentified specimens; NCBI: missing
`lat_lon`).

## Usage

``` r
extract_edna_occurrences(ncbi_results, bold_raw)
```

## Arguments

- ncbi_results:

  Tibble as returned by
  [`fetch_ncbi_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_ncbi_sequences.md)'s
  `results` (must have `organism`, `lat_lon`, `collection_date`,
  `accession` columns)

- bold_raw:

  Character scalar, raw BOLD TSV as returned by
  [`fetch_bold_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_bold_sequences.md)

## Value

Tibble with columns `source` (`"NCBI"`/`"BOLD"`), `species`, `lon`,
`lat`, `date` (`Date`, `NA` where uncollectable), `record_id`
