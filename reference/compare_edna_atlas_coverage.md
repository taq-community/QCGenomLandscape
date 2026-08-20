# Compare eDNA vs traditional occurrence coverage against the Atlas, on a hex grid

The heavy lifting: classifies each Atlas record as `"eDNA"` (one of
[`atlas_molecular_datasets()`](https://taq-community.github.io/QCGenomLandscape/reference/atlas_molecular_datasets.md))
or `"traditional"`, spatially bins every record into `hex_grid` via a
DuckDB spatial join (fast even at Atlas scale – tens of millions of rows
join against a ~18k-cell grid in well under a minute), and folds in
`edna_occurrences` (e.g. from
[`extract_edna_occurrences()`](https://taq-community.github.io/QCGenomLandscape/reference/extract_edna_occurrences.md))
as additional eDNA evidence, binned into the same grid via
[`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html).

## Usage

``` r
compare_edna_atlas_coverage(
  atlas_parquet_path,
  hex_grid,
  edna_occurrences,
  bdqc_species = NULL,
  molecular_datasets = atlas_molecular_datasets(),
  con = NULL
)
```

## Arguments

- atlas_parquet_path:

  Character, local path to an Atlas parquet export, e.g. from
  [`fetch_atlas_parquet()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_atlas_parquet.md)

- hex_grid:

  `sf` polygon grid with a `cell_id` column, e.g. from
  [`build_hex_grid_cells()`](https://taq-community.github.io/QCGenomLandscape/reference/build_hex_grid_cells.md)

- edna_occurrences:

  Tibble with columns `species`, `lon`, `lat`, `date`, e.g. from
  [`extract_edna_occurrences()`](https://taq-community.github.io/QCGenomLandscape/reference/extract_edna_occurrences.md)

- bdqc_species:

  Character vector of species to restrict the comparison to, or `NULL`
  for no restriction (default `NULL`)

- molecular_datasets:

  Character vector of Atlas `dataset_name` values to treat as eDNA,
  default
  [`atlas_molecular_datasets()`](https://taq-community.github.io/QCGenomLandscape/reference/atlas_molecular_datasets.md)

- con:

  An existing DBI connection to a DuckDB database, or `NULL` (default)
  to open and close a temporary one; injectable for testing

## Value

A list with `edna_agg` and `traditional_agg` tibbles, each with columns
`species`, `cell_id`, `year_obs`, `n`

## Details

No deduplication is attempted between `edna_occurrences` and Atlas's own
molecular datasets (some of which, e.g. iBOL, share underlying records
with a BOLD-sourced `edna_occurrences`) – this has little effect on the
presence-based metrics
[`summarize_edna_contribution()`](https://taq-community.github.io/QCGenomLandscape/reference/summarize_edna_contribution.md)
computes from the result, which count distinct (species, cell) /
(species, year) pairs, not raw record volume.
