# Drop subspecies/infraspecific records from NCBI results

NCBI's `[Organism]` search matches a species name *and* everything below
it taxonomically (subspecies, breeds, domestic forms), e.g. querying
`"Canis lupus[Organism]"` also returns `"Canis lupus familiaris"`
records. Left unfiltered, these can dominate a species' sequence count
(domestic dog barcodes outnumber wild wolf ones ~3:1 in the raw NCBI
pull) and bias any coverage/quality analysis keyed on the BDQC species
list. Keeping only rows whose `organism` exactly matches a name in
`species` removes them.

## Usage

``` r
filter_named_species(results, species)
```

## Arguments

- results:

  Tibble with an `organism` column, e.g.
  `fetch_ncbi_sequences()$results`

- species:

  Character vector of exact species names to keep (typically the BDQC
  species list)

## Value

`results`, filtered to rows where `organism %in% species`
