# Flag whether a GenBank record's DEFINITION describes a complete genome

GenBank distinguishes `"complete genome"` from `"partial genome"` in its
own DEFINITION wording (e.g.
`"Boreogadus saida ... mitochondrion, complete genome."` vs
`"... mitochondrion, partial genome."`) – the authoritative signal for
"is this a full genome", direct from the submitter, rather than
inferring it from how many `/gene=` tags happen to be annotated on the
record (a multi-gene *fragment* spanning just two adjacent genes would
also have multiple tags without being anywhere near a full genome).

## Usage

``` r
is_complete_genome(definition)
```

## Arguments

- definition:

  Character vector, e.g.
  [`parse_gb_records()`](https://taq-community.github.io/QCGenomLandscape/reference/parse_gb_records.md)'s
  `definition` column

## Value

Logical vector, same length as `definition`; `NA` where `definition` is
`NA`

## Examples

``` r
is_complete_genome(c(
  "Boreogadus saida mitochondrion, complete genome.",
  "Boreogadus saida mitochondrion, partial genome.",
  "Boreogadus saida cytochrome oxidase subunit I (COI) gene, partial cds."
))
#> [1]  TRUE FALSE FALSE
```
