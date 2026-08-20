# Parse a BOLD-style `coord` string into decimal degrees

Parse a BOLD-style `coord` string into decimal degrees

## Usage

``` r
parse_bold_coord(coord_string)
```

## Arguments

- coord_string:

  A string like `"[45.45, -75.77]"` (`"[lat, lon]"`)

## Value

A list with `lat` and `lon` as numeric values (`NA_real_` if the string
is missing, empty, or doesn't contain exactly two numbers)

## Examples

``` r
parse_bold_coord("[45.45, -75.77]")
#> $lat
#> [1] 45.45
#> 
#> $lon
#> [1] -75.77
#> 
parse_bold_coord(NA)
#> $lat
#> [1] NA
#> 
#> $lon
#> [1] NA
#> 
```
