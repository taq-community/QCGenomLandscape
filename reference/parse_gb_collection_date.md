# Parse a GenBank/BOLD-style collection-date string into a `Date`

Collection dates in GenBank (and, less often, BOLD) records come in
several formats of varying precision – `"DD-Mon-YYYY"` (the common
case), `"Mon-YYYY"`/`"YYYY-MM"`/`"YYYY"` (partial dates, normalized to
the first of the month/year), plus ISO `"YYYY-MM-DD"`. Occasionally the
field holds non-date junk instead (e.g. a locality string or a
coordinate pair that leaked into it upstream); that returns `NA`.

## Usage

``` r
parse_gb_collection_date(x)
```

## Arguments

- x:

  Character vector of collection-date values

## Value

`Date` vector, the same length as `x`

## Details

[`as.Date()`](https://rdrr.io/r/base/as.Date.html) is not used directly
here: with an explicit `format=` it does *partial*, not full-string,
matching (`as.Date("58.298 N 8.539 E", format = "%Y")` silently parses
to the year 58 using today's month/day) and formats lacking `%d`
(`"%b-%Y"`, `"%Y-%m"`) fail outright rather than defaulting the day.
Each recognized pattern is instead converted to an explicit ISO string
first, so the final parse is unambiguous.

## Examples

``` r
parse_gb_collection_date(c("28-Sep-2014", "Sep-2023", "2020", "junk", NA))
#> [1] "2014-09-28" "2023-09-01" "2020-01-01" NA           NA          
```
