# Fetch BOLD Systems records for a query, optionally saving to disk

Fetch BOLD Systems records for a query, optionally saving to disk

## Usage

``` r
fetch_bold_sequences(
  query = "geo:province/state:Quebec",
  extent = "full",
  out_path = NULL,
  request_fn = httr2::request
)
```

## Arguments

- query:

  Character, BOLD query string, default `"geo:province/state:Quebec"`

- extent:

  Character, default `"full"`

- out_path:

  Character or `NULL`; if not `NULL`, writes the raw TSV there

- request_fn:

  Function with signature `(url)` returning an `httr2` request, default
  [`httr2::request()`](https://httr2.r-lib.org/reference/request.html);
  injectable for testing without network access

## Value

Invisibly, the raw TSV text (character scalar)
