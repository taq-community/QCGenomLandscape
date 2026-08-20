# Summarize an eDNA-vs-traditional coverage comparison into contribution metrics

Pure post-processing over the output of
[`compare_edna_atlas_coverage()`](https://taq-community.github.io/QCGenomLandscape/reference/compare_edna_atlas_coverage.md)
– no I/O, so it's cheap to recompute (e.g. per taxonomic subgroup)
without re-running the spatial join.

## Usage

``` r
summarize_edna_contribution(edna_agg, traditional_agg)
```

## Arguments

- edna_agg, traditional_agg:

  Tibbles with columns `species`, `cell_id`, `year_obs`, `n`, e.g. from
  [`compare_edna_atlas_coverage()`](https://taq-community.github.io/QCGenomLandscape/reference/compare_edna_atlas_coverage.md)

## Value

A list: `edna_species`, `traditional_species`, `edna_only_species`,
`traditional_only_species`, `shared_species` (character vectors);
`novel_cells`, `novel_years` (tibbles of `species` + the novel
cell/year); and summary counts/percentages (`n_species_edna`,
`n_species_edna_only`, `pct_species_edna_only`, `n_species_traditional`,
`n_species_traditional_only`, `pct_species_traditional_only`,
`n_species_shared`, `jaccard_species` (0 = disjoint species sets, 1 =
identical), `n_cells_edna`, `n_cells_novel`, `pct_cells_novel`,
`n_years_edna`, `n_years_novel`, `pct_years_novel`)
