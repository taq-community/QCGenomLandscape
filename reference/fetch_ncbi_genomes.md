# Query full-genome availability for a vector of species, batched

Splits `species` into batches of `batch_size` and resolves nuclear and
mitochondrial genome availability independently per batch – see
`resolve_nuclear_genomes()` and `resolve_mitochondrial_genomes()` for
why they use different strategies (nuclear hits are rare enough for a
fast existence-check path to pay off; mitochondrial hits are not, so
that dimension goes straight to a batched search + summary attribution).

## Usage

``` r
fetch_ncbi_genomes(
  species,
  batch_size = 25,
  search_fn = rentrez::entrez_search,
  summary_fn = rentrez::entrez_summary,
  progress = TRUE
)
```

## Arguments

- species:

  Character vector of species names (deduplicated internally)

- batch_size:

  Integer, species OR'd together per query, default 25. Keep well under
  NCBI's query-length limits.

- search_fn:

  Function with signature `(db, term, retmax)`, default
  [`rentrez::entrez_search()`](https://docs.ropensci.org/rentrez/reference/entrez_search.html);
  injectable for testing without network access

- summary_fn:

  Function with signature `(db, id)`, default
  [`rentrez::entrez_summary()`](https://docs.ropensci.org/rentrez/reference/entrez_summary.html)

- progress:

  Logical, show a progress bar, default `TRUE`

## Value

Tibble, one row per species: `species`, `nuclear_genome`,
`mitochondrial_genome`, `nuclear_count`, `mitochondrial_count`,
`nuclear_accessions`, `mitochondrial_accessions`, `error`
