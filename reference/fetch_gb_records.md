# Fetch raw GenBank flat-file text for a set of accessions

Fetch raw GenBank flat-file text for a set of accessions

## Usage

``` r
fetch_gb_records(
  accessions,
  batch_size = 50,
  fetch_fn = rentrez::entrez_fetch,
  progress = TRUE
)
```

## Arguments

- accessions:

  Character vector of NCBI accession numbers

- batch_size:

  Integer, number of accessions fetched per request, default 50

- fetch_fn:

  Function with signature `(db, id, rettype, retmode)`, default
  [`rentrez::entrez_fetch()`](https://docs.ropensci.org/rentrez/reference/entrez_fetch.html);
  injectable for testing without network access

- progress:

  Logical, show a progress bar, default `TRUE`

## Value

A list of character scalars, one raw GenBank flat-file blob per batch
