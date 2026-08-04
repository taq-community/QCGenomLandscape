# Parse accession/gene/sequence out of GenBank flat-file text

Parse accession/gene/sequence out of GenBank flat-file text

## Usage

``` r
parse_gb_records(gb_text)
```

## Arguments

- gb_text:

  Character scalar containing one or more `"//"`-delimited GenBank
  flat-file records

## Value

Tibble with columns `accession`, `gene` (`;`-separated if a record has
multiple `/gene=` tags), `sequence` (lowercase-normalized `acgt` only,
`NA` if the record has no `ORIGIN` block)
