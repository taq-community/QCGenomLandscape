# Parse an NCBI-style lat/lon string into decimal degrees

Parse an NCBI-style lat/lon string into decimal degrees

## Usage

``` r
parse_latlon(latlon_string)
```

## Arguments

- latlon_string:

  A string like "45.5 N 73.6 W" or "45.5N 73.6W"

## Value

A list with `lat` and `lon` as numeric values (`NA_real_` if the string
is missing, empty, or doesn't contain at least two numbers)

## Examples

``` r
parse_latlon("45.5 N 73.6 W")
#> $lat
#> [1] 45.5
#> 
#> $lon
#> [1] -73.6
#> 
parse_latlon("33.9 S 18.4 E")
#> $lat
#> [1] -33.9
#> 
#> $lon
#> [1] 18.4
#> 
parse_latlon(NA)
#> $lat
#> [1] NA
#> 
#> $lon
#> [1] NA
#> 
```
