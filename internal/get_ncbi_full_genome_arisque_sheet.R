#!/usr/bin/env Rscript
# Query NCBI full-genome availability for species listed in the private
# "À risque" Google Sheet, instead of the full BDQC species list (that's
# what the main _targets.R pipeline's ncbi_genomes target does). Reuses the
# same fetch_ncbi_genomes() the package exports -- only the species source
# differs.
library(QCGenomLandscape)

SHEET_ID <- "1W5__6vlrnY1IeZPYEeTAx7a8bG8oLTpj4GRISthXsos"

log_dir <- "logs"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

log_file <- file.path(log_dir, sprintf("ncbi_full_genome_arisque_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
logger::log_appender(logger::appender_file(log_file))
logger::log_threshold(logger::INFO)

logger::log_info("Logger initialized. Writing to: {log_file}")
rentrez::set_entrez_key(Sys.getenv("NCBI_API_KEY"))

# Species at risk from Google Sheet -- "À risque" tab
googlesheets4::gs4_auth()
qc_species <- googlesheets4::read_sheet(SHEET_ID, sheet = "À risque") |>
  dplyr::pull(`Espèce`) |>
  unique() |>
  na.omit()

logger::log_info("Starting genome queries for {length(qc_species)} species at risk")

genome_results <- fetch_ncbi_genomes(qc_species)

saveRDS(genome_results, "results/ncbi_full_genome_species_at_risk.rds")

logger::log_success("Completed!")
