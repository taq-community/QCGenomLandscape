# Fetch raw GenBank flat-file text for a set of accessions

Transient network/SSL errors (e.g. NCBI dropping the connection mid-
response) are retried with backoff; a batch that still fails after
`max_retries` is logged and dropped rather than failing the whole call –
with batch sizes in the hundreds, one flaky batch shouldn't lose every
other batch's sequences.

## Usage

``` r
fetch_gb_records(
  accessions,
  batch_size = 50,
  fetch_fn = rentrez::entrez_fetch,
  max_retries = 3,
  sleep_fn = Sys.sleep,
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

- max_retries:

  Integer, attempts per batch before giving up on it, default 3
  (exponential backoff between attempts: 2s, 4s, ...)

- sleep_fn:

  Function with signature `(seconds)`, default
  [`Sys.sleep()`](https://rdrr.io/r/base/Sys.sleep.html); injectable so
  retry backoff doesn't slow down tests

- progress:

  Logical, show a progress bar, default `TRUE`

## Value

A list of character scalars, one raw GenBank flat-file blob per
successful batch (shorter than the number of batches if any were dropped
after exhausting retries)
