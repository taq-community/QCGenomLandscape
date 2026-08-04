# internal/

One-off reporting scripts tied to a specific private Google Sheet (the "À risque" species-at-risk tracker), not part of the reproducible `_targets.R` pipeline and not exported package functions.

- `get_ncbi_full_genome_arisque_sheet.R` — same genome-availability query as the main pipeline, but sourcing the species list from the Sheet instead of the full BDQC list.
- `fill_arisque_sheet.R` — builds the species-at-risk summary table (COSEPAC/LEMV status, qPCR kit availability, genome-sequencing status) and writes it back to the Sheet's "À risque" tab.

Requirements not declared as package dependencies (install separately): `googlesheets4`, `taxadb` (needs `taxadb::td_create("itis")` run once locally). Both scripts need Google Sheets access to the sheet ID hardcoded in each file, and `googlesheets4::gs4_auth()` will prompt for OAuth on first run.
