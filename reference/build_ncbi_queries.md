# Build NCBI search query strings from a BDQC species list and primer map

With `batch_size > 1`, species sharing the same `query_marker` are
combined into a single OR'd query
(`"(Sp1[Organism] OR Sp2[Organism] OR ...) AND marker[Gene]"`), cutting
the number of Entrez searches by roughly `batch_size`x. Species are no
longer recoverable by parsing the query string in that case – use the
`organism` field NCBI returns in
[`fetch_ncbi_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_ncbi_sequences.md)'s
results instead (already present in its output).

## Usage

``` r
build_ncbi_queries(species_df, query_primers, voucher = TRUE, batch_size = 1)
```

## Arguments

- species_df:

  Data frame from the BDQC species list csv, already filtered to
  `rank == "species"` (must have `species` and `group_en` columns)

- query_primers:

  Data frame from the primers-map csv (must have `group` and
  `query_marker` columns)

- voucher:

  Logical; if `TRUE` (default), appends `AND voucher[Title]` to each
  query so results are restricted to specimen-voucher-backed records

- batch_size:

  Integer, how many species (sharing the same marker) to OR together per
  query, default 1 (one query per species, matching the original
  per-species scripts). Keep this well under NCBI's query-length limits
  – 20-50 is a reasonable range in practice.

## Value

Character vector of unique Entrez query strings
