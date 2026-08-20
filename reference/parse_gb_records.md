# Parse accession/definition/gene/sequence out of GenBank flat-file text

Parse accession/definition/gene/sequence out of GenBank flat-file text

## Usage

``` r
parse_gb_records(gb_text)
```

## Arguments

- gb_text:

  Character scalar containing one or more `"//"`-delimited GenBank
  flat-file records

## Value

Tibble with columns `accession`, `definition` (the GenBank `DEFINITION`
field, whitespace-collapsed to one line – see
[`is_complete_genome()`](https://taq-community.github.io/QCGenomLandscape/reference/is_complete_genome.md)),
`gene` (`;`-separated if a record has multiple `/gene=` tags),
`sequence` (lowercase-normalized `acgt` only, `NA` if the record has no
`ORIGIN` block)
