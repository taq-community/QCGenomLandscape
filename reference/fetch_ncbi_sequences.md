# Fetch NCBI nucleotide records for a set of species/marker queries

The batching/error-handling loop shared by the voucher and non-voucher
NCBI query scripts. `search_fn`/`summary_fn` are injectable so the
batching, high-ID-count flagging, and error-accumulation logic can be
unit tested without live network access or an `NCBI_API_KEY`.

## Usage

``` r
fetch_ncbi_sequences(
  queries,
  retmax = 5000,
  batch_size = 200,
  high_id_threshold = 500,
  search_fn = rentrez::entrez_search,
  summary_fn = rentrez::entrez_summary,
  progress = TRUE
)
```

## Arguments

- queries:

  Character vector of Entrez query strings, e.g. from
  [`build_ncbi_queries()`](https://taq-community.github.io/QCGenomLandscape/reference/build_ncbi_queries.md)

- retmax:

  Integer, per-search max IDs to retrieve, default 5000

- batch_size:

  Integer, `entrez_summary` batch size, default 200

- high_id_threshold:

  Integer, queries returning more IDs than this are flagged in
  `high_id_queries`, default 500

- search_fn:

  Function with signature `(db, term, retmax)`, default
  [`rentrez::entrez_search()`](https://docs.ropensci.org/rentrez/reference/entrez_search.html)

- summary_fn:

  Function with signature `(db, id)`, default
  [`rentrez::entrez_summary()`](https://docs.ropensci.org/rentrez/reference/entrez_summary.html)

- progress:

  Logical, show a progress bar, default `TRUE`

## Value

A list with three elements: `results` (tibble of parsed sequence
summaries), `deficient_queries` (list of queries that errored), and
`high_id_queries` (list of queries whose ID count exceeded
`high_id_threshold`)
