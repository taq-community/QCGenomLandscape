# Stop the pipeline when Entrez queries were lost

Called as its own `targets` target so that an incomplete NCBI corpus is
a build failure with the affected species named, rather than a silent
hole that only shows up as an implausible zero in a figure. Set
`allow_deficient = TRUE` (or `QCGENOM_ALLOW_DEFICIENT=true`) to
downgrade it to a warning when a partial corpus is genuinely acceptable.

## Usage

``` r
assert_no_deficient_queries(
  deficient_queries,
  species,
  allow_deficient = identical(tolower(Sys.getenv("QCGENOM_ALLOW_DEFICIENT")), "true")
)
```

## Arguments

- deficient_queries:

  List as returned in `fetch_ncbi_sequences()$deficient_queries`

- species:

  Character vector of species names, typically the BDQC list

- allow_deficient:

  Logical, warn instead of erroring

## Value

The tibble from
[`deficient_query_species()`](https://taq-community.github.io/QCGenomLandscape/reference/deficient_query_species.md),
invisibly
