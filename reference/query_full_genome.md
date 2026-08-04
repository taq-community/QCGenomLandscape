# Query NCBI genome/nucleotide databases for full-genome availability

Query NCBI genome/nucleotide databases for full-genome availability

## Usage

``` r
query_full_genome(
  species_name,
  query_index = 1L,
  total_queries = 1L,
  search_fn = rentrez::entrez_search
)
```

## Arguments

- species_name:

  Character scalar, species name

- query_index:

  Integer, for logging only, default 1

- total_queries:

  Integer, for logging only, default 1

- search_fn:

  Function with signature `(db, term, retmax)`, default
  [`rentrez::entrez_search()`](https://docs.ropensci.org/rentrez/reference/entrez_search.html);
  injectable for testing without network access

## Value

Tibble with one row: `species`, `nuclear_genome`,
`mitochondrial_genome`, `nuclear_count`, `mitochondrial_count`,
`nuclear_accessions`, `mitochondrial_accessions`, `error`
