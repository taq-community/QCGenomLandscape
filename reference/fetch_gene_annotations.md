# Fetch gene annotations for a set of accessions via GenBank XML

Fetch gene annotations for a set of accessions via GenBank XML

## Usage

``` r
fetch_gene_annotations(
  accessions,
  batch_size = 200,
  fetch_fn = rentrez::entrez_fetch,
  progress = TRUE
)
```

## Arguments

- accessions:

  Character vector of NCBI accession numbers

- batch_size:

  Integer, number of accessions fetched per request, default 200

- fetch_fn:

  Function with signature `(db, id, rettype, retmode)`, default
  [`rentrez::entrez_fetch()`](https://docs.ropensci.org/rentrez/reference/entrez_fetch.html);
  injectable for testing without network access

- progress:

  Logical, show a progress bar, default `TRUE`

## Value

Tibble with columns `accession`, `gene`, `location`
