# Load and normalize a conservation-risk-status table

Unifies the CA/COSEPAC and QC/LEMV risk-status loading blocks previously
duplicated between `create_dataframe.R` and `taxon_representation.R`.

## Usage

``` r
load_risk_status(path, jurisdiction = c("CA", "QC"), translate = FALSE)
```

## Arguments

- path:

  Character, CSV path (`data/CA_especes_en_peril.csv` for `"CA"`,
  `data/QC_especes_en_peril.csv` for `"QC"`)

- jurisdiction:

  `"CA"` or `"QC"`

- translate:

  Logical; if `TRUE`, recodes French status labels to English (as
  `taxon_representation.R` did); default `FALSE` keeps the original
  French labels (as `create_dataframe.R` did)

## Value

Tibble with columns `species`, `status`, `jurisdiction`
