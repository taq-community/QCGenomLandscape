# Retry an Entrez call across transient failures

Entrez fails transiently often enough that a single attempt is not a
reliable read: DNS resolution times out, the load balancer returns 502
on large batches, and rentrez occasionally throws
`"subscript out of bounds"` parsing a truncated response. Every one of
those is retryable, and every one of them silently cost this pipeline a
whole batch of species before this helper existed – the 2026-08-19 run
lost all 24 Quebec amphibians to a single 10-second DNS timeout on query
1 of 954.

## Usage

``` r
with_entrez_retry(
  expr,
  max_attempts = 4,
  sleep_fn = Sys.sleep,
  label = "entrez call"
)
```

## Arguments

- expr:

  Quoted expression performing the Entrez call

- max_attempts:

  Integer, total attempts before giving up, default 4

- sleep_fn:

  Function with signature `(seconds)`, default
  [`Sys.sleep()`](https://rdrr.io/r/base/Sys.sleep.html); injectable so
  retry backoff doesn't slow down tests

- label:

  Character, short description used in log messages

## Value

The value of `expr`, or a `condition` object if every attempt failed
