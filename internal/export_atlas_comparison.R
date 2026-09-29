#!/usr/bin/env Rscript
# Promote the eDNA-vs-Atlas comparison out of the knitr cache into results/.
#
# compare_edna_atlas_coverage() runs against the 2.1 GB Atlas parquet and
# lives in a `cache: true` chunk of portrait-genomique-bdqc.qmd, so the only
# copy of its output was a knitr cache keyed by a chunk hash -- not something
# the manuscript can depend on. This writes the same objects to a stable file
# so docs/manuscript.qmd can compute Table 1 and Figures 1-3 from them
# instead of carrying transcribed numbers.
#
# Usage: Rscript internal/export_atlas_comparison.R

suppressMessages(library(logger))

cache_dir <- "portrait-genomique-bdqc_cache/html"
db <- list.files(cache_dir, pattern = "^atlas-compare_.*\\.rdb$", full.names = TRUE)
if (length(db) != 1) {
  stop(sprintf("expected exactly one atlas-compare cache db, found %d", length(db)))
}

e <- new.env()
lazyLoad(sub("\\.rdb$", "", db), envir = e)
log_info("Loaded from cache: {paste(ls(e), collapse = ', ')}")

out <- list(
  contribution = e$contribution,
  coverage = e$coverage,
  edna_occ = e$edna_occ,
  species_groupe = e$species_groupe,
  hex_grid = e$hex_grid,
  qc_boundary = e$qc_boundary,
  cache_file = basename(db),
  exported_at = Sys.time()
)

dir.create("results", showWarnings = FALSE)
saveRDS(out, "results/edna_atlas_comparison.rds")
log_success(
  "results/edna_atlas_comparison.rds -- {out$contribution$n_species_edna} eDNA species, ",
  "{nrow(out$coverage$edna_agg)} eDNA (species, cell, year) rows"
)
