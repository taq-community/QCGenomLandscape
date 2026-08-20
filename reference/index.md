# Package index

## All functions

- [`align_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/align_sequences.md)
  : Align a set of DNA sequences

- [`assign_gene_group()`](https://taq-community.github.io/QCGenomLandscape/reference/assign_gene_group.md)
  : Classify a gene name into a coarse marker group

- [`atlas_molecular_datasets()`](https://taq-community.github.io/QCGenomLandscape/reference/atlas_molecular_datasets.md)
  : The 16 Atlas datasets whose records are already
  molecular/eDNA-derived

- [`build_hex_grid_cells()`](https://taq-community.github.io/QCGenomLandscape/reference/build_hex_grid_cells.md)
  : Build a full hexagonal grid over a boundary, with stable cell IDs

- [`build_ncbi_queries()`](https://taq-community.github.io/QCGenomLandscape/reference/build_ncbi_queries.md)
  : Build NCBI search query strings from a BDQC species list and primer
  map

- [`build_sequence_qc_table()`](https://taq-community.github.io/QCGenomLandscape/reference/build_sequence_qc_table.md)
  : Build the per-record sequence-QC table from fetched GenBank batches

- [`build_summary_dataframe()`](https://taq-community.github.io/QCGenomLandscape/reference/build_summary_dataframe.md)
  : Build the per-species genomic-data summary table

- [`build_summary_dataframe_from_files()`](https://taq-community.github.io/QCGenomLandscape/reference/build_summary_dataframe_from_files.md)
  : Build the per-species genomic-data summary table, reading inputs
  from disk

- [`classify_taxon_group()`](https://taq-community.github.io/QCGenomLandscape/reference/classify_taxon_group.md)
  : Classify BDQC taxonomic groups into finer display categories

- [`compare_edna_atlas_coverage()`](https://taq-community.github.io/QCGenomLandscape/reference/compare_edna_atlas_coverage.md)
  : Compare eDNA vs traditional occurrence coverage against the Atlas,
  on a hex grid

- [`extract_edna_occurrences()`](https://taq-community.github.io/QCGenomLandscape/reference/extract_edna_occurrences.md)
  : Normalize NCBI + BOLD records into a common eDNA occurrence table

- [`fetch_atlas_parquet()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_atlas_parquet.md)
  : Fetch a Biodiversité Québec Atlas public-data export, caching it
  locally

- [`fetch_bold_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_bold_sequences.md)
  : Fetch BOLD Systems records for a query, optionally saving to disk

- [`fetch_gb_records()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_gb_records.md)
  : Fetch raw GenBank flat-file text for a set of accessions

- [`fetch_gene_annotations()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_gene_annotations.md)
  : Fetch gene annotations for a set of accessions via GenBank XML

- [`fetch_ncbi_genomes()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_ncbi_genomes.md)
  : Query full-genome availability for a vector of species, batched

- [`fetch_ncbi_sequences()`](https://taq-community.github.io/QCGenomLandscape/reference/fetch_ncbi_sequences.md)
  : Fetch NCBI nucleotide records for a set of species/marker queries

- [`flag_barcode_gap_outliers()`](https://taq-community.github.io/QCGenomLandscape/reference/flag_barcode_gap_outliers.md)
  : Flag DNA sequences whose distance to conspecifics is a group outlier

- [`flag_length_outliers()`](https://taq-community.github.io/QCGenomLandscape/reference/flag_length_outliers.md)
  : Flag sequence-length outliers within each marker/gene group

- [`flag_within_boundary()`](https://taq-community.github.io/QCGenomLandscape/reference/flag_within_boundary.md)
  : Flag points as within a polygon boundary

- [`has_stop_codon_coi()`](https://taq-community.github.io/QCGenomLandscape/reference/has_stop_codon_coi.md)
  : Flag DNA sequences with no clean reading frame (COI / SGC4 code)

- [`is_complete_genome()`](https://taq-community.github.io/QCGenomLandscape/reference/is_complete_genome.md)
  : Flag whether a GenBank record's DEFINITION describes a complete
  genome

- [`load_canvec_boundary()`](https://taq-community.github.io/QCGenomLandscape/reference/load_canvec_boundary.md)
  : Load a CanVec administrative-boundary layer, optionally filtered

- [`load_risk_status()`](https://taq-community.github.io/QCGenomLandscape/reference/load_risk_status.md)
  : Load and normalize a conservation-risk-status table

- [`make_hex_grid()`](https://taq-community.github.io/QCGenomLandscape/reference/make_hex_grid.md)
  : Build a hexagonal summary grid of point counts within a boundary

- [`parse_bold_coord()`](https://taq-community.github.io/QCGenomLandscape/reference/parse_bold_coord.md)
  :

  Parse a BOLD-style `coord` string into decimal degrees

- [`parse_gb_collection_date()`](https://taq-community.github.io/QCGenomLandscape/reference/parse_gb_collection_date.md)
  :

  Parse a GenBank/BOLD-style collection-date string into a `Date`

- [`parse_gb_records()`](https://taq-community.github.io/QCGenomLandscape/reference/parse_gb_records.md)
  : Parse accession/definition/gene/sequence out of GenBank flat-file
  text

- [`parse_latlon()`](https://taq-community.github.io/QCGenomLandscape/reference/parse_latlon.md)
  : Parse an NCBI-style lat/lon string into decimal degrees

- [`plot_gene_prevalence()`](https://taq-community.github.io/QCGenomLandscape/reference/plot_gene_prevalence.md)
  : Donut chart: gene prevalence by taxonomic group

- [`plot_risk_status_coverage()`](https://taq-community.github.io/QCGenomLandscape/reference/plot_risk_status_coverage.md)
  : Donut chart: genomic-data coverage by species-at-risk status

- [`query_full_genome()`](https://taq-community.github.io/QCGenomLandscape/reference/query_full_genome.md)
  : Query NCBI genome/nucleotide databases for full-genome availability

- [`score_sequence_quality()`](https://taq-community.github.io/QCGenomLandscape/reference/score_sequence_quality.md)
  : Compute basic sequence-quality metrics

- [`summarize_edna_contribution()`](https://taq-community.github.io/QCGenomLandscape/reference/summarize_edna_contribution.md)
  : Summarize an eDNA-vs-traditional coverage comparison into
  contribution metrics

- [`swift_auth()`](https://taq-community.github.io/QCGenomLandscape/reference/swift_auth.md)
  : Authenticate against an OpenStack Keystone endpoint

- [`upload_to_swift()`](https://taq-community.github.io/QCGenomLandscape/reference/upload_to_swift.md)
  : Upload a local file to Arbutus (OpenStack Swift) object storage
