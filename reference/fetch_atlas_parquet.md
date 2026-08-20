# Fetch a Biodiversité Québec Atlas public-data export, caching it locally

Downloads via the same public, tokenless S3 bucket documented at
`source("http://atlas.biodiversite-quebec.ca/bq-atlas-parquet.R")` –
re-implemented directly (rather than
[`source()`](https://rdrr.io/r/base/source.html)-ing that URL) so the
fetch is inspectable and testable. Files are ~2 GB; already-downloaded
exports are reused, not re-fetched.

## Usage

``` r
fetch_atlas_parquet(
  export_date = NULL,
  dest_dir = "data/atlas",
  base_url = "https://object-arbutus.cloud.computecanada.ca/bq-io/atlas/parquet/",
  list_dates_fn = function(url) {
     dates <- utils::read.csv(url, header = FALSE,
    col.names = "date")$date
dates[nzchar(trimws(dates))]
 },
  download_fn = function(url, destfile) {
     utils::download.file(url, destfile, mode =
    "wb", quiet = FALSE)
 }
)
```

## Arguments

- export_date:

  Character `"YYYY-MM-DD"`, default the latest export listed in the
  bucket's `atlas_export_dates.csv`

- dest_dir:

  Directory to cache the downloaded parquet file in, default
  `"data/atlas"`

- base_url:

  Base URL of the Atlas parquet exports

- list_dates_fn:

  Function with signature `(url)` returning the export dates as a
  character vector, default reads `atlas_export_dates.csv`; injectable
  for testing without network access

- download_fn:

  Function with signature `(url, destfile)`, default
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html);
  injectable for testing without network access

## Value

Character scalar, local path to the parquet file
