# Build the per-species genomic-data summary table, reading inputs from disk

Thin I/O wrapper around
[`build_summary_dataframe()`](https://taq-community.github.io/QCGenomLandscape/reference/build_summary_dataframe.md)
for interactive/ad hoc use outside the `_targets.R` pipeline (which
calls
[`build_summary_dataframe()`](https://taq-community.github.io/QCGenomLandscape/reference/build_summary_dataframe.md)
directly with in-memory pipeline objects).

## Usage

``` r
build_summary_dataframe_from_files(
  ncbi_path = "results/ncbi_results.rds",
  genes_path = "results/genes_subsamp_50_df.rds",
  bdqc_path = "data/bdqc_list_01122025.csv",
  ca_risk_path = "data/CA_especes_en_peril.csv",
  qc_risk_path = "data/QC_especes_en_peril.csv"
)
```

## Arguments

- ncbi_path:

  Character, path to `fetch_ncbi_sequences()$results`, saved as RDS

- genes_path:

  Character, path to
  [`fetch_gene_annotations()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_gene_annotations.md)
  output, saved as RDS

- bdqc_path:

  Character, path to the BDQC species list CSV

- ca_risk_path, qc_risk_path:

  Character, paths to the CA/QC risk-status CSVs

## Value

Tibble, one row per species (see
[`build_summary_dataframe()`](https://taq-community.github.io/QCGenomLandscape/reference/build_summary_dataframe.md))
