#!/usr/bin/env Rscript
# Recover the Entrez batches lost to transient errors in the 2026-08-19 run.
#
# Five of 954 queries failed and were recorded in
# `ncbi_sequences$deficient_queries` without anything downstream reading that
# record. Query 1 -- a DNS timeout -- took all 24 Quebec amphibians with it and
# was published as "Amphibians: 0% coverage" in Figure 4. Three lichen/fungal
# batches hit "subscript out of bounds" and one hit a 502.
#
# This re-runs exactly those five queries through the now-retrying
# fetch_ncbi_sequences(), merges the recovered records into
# results/ncbi_results.rds, and refreshes the gene annotations and sequence-QC
# rows for the newly covered species. It is a one-shot repair, not part of the
# pipeline: a full `tar_make()` with the fixed code would produce the same
# corpus from scratch.
#
# Usage: Rscript internal/refetch_deficient_queries.R

suppressMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(logger)
})

for (f in c("R/ncbi_utils.R", "R/ncbi_sequences.R", "R/genbank.R",
            "R/gene_classification.R", "R/utils.R")) {
  if (file.exists(f)) source(f)
}

stopifnot(dir.exists("results"))
ncbi_results <- readRDS("results/ncbi_results.rds")
deficient <- readRDS("_targets/objects/ncbi_sequences")$deficient_queries
bdqc <- read.csv("data/bdqc_list_01122025.csv") |>
  dplyr::filter(rank == "species")

queries <- vapply(deficient, function(d) d$query, character(1))
log_info("Re-running {length(queries)} failed queries")

recovered <- fetch_ncbi_sequences(queries, progress = FALSE)
if (length(recovered$deficient_queries)) {
  log_error("{length(recovered$deficient_queries)} queries still failing")
}
new_records <- filter_named_species(recovered$results, bdqc$species)
log_info("Recovered {nrow(new_records)} records for {dplyr::n_distinct(new_records$organism)} species")

merged <- dplyr::bind_rows(ncbi_results, new_records) |>
  dplyr::distinct(uid, .keep_all = TRUE)
log_info(
  "Corpus: {nrow(ncbi_results)} -> {nrow(merged)} records, ",
  "{dplyr::n_distinct(ncbi_results$organism)} -> {dplyr::n_distinct(merged$organism)} species"
)

saveRDS(merged, "results/ncbi_results.rds")
saveRDS(new_records, "results/recovered_records.rds")
log_success("Wrote results/ncbi_results.rds")
