#!/usr/bin/env Rscript
# Bring the derived results/ objects back in line with the repaired corpus.
#
# Companion to internal/refetch_deficient_queries.R. Two independent repairs
# land here:
#
#   1. the 5,680 records recovered from the five failed Entrez batches need
#      sequence-QC rows and gene annotations of their own;
#   2. `gene_annotations` predates the DEFINITION fallback and carries no
#      `definition` column, so every ITS record in it is still unclassified.
#      The definitions already exist in sequence_qc.rds (parse_gb_records
#      captures them for every accession), so they are joined in rather than
#      re-fetched.
#
# Usage: Rscript internal/refresh_after_recovery.R

suppressMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(logger)
})

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

strip_version <- function(x) sub("\\.\\d+$", "", x)

ncbi <- readRDS("results/ncbi_results.rds")
seq_data <- readRDS("results/sequence_qc.rds")
# genes_saved writes results/genes_subsamp_50_df.rds, but the 2026-08-19 run
# left only the targets copy on disk -- read whichever exists.
genes_path <- "results/genes_subsamp_50_df.rds"
genes <- if (file.exists(genes_path)) {
  readRDS(genes_path)
} else {
  readRDS("_targets/objects/gene_annotations")
}
bdqc <- read.csv("data/bdqc_list_01122025.csv") |> dplyr::filter(rank == "species")

# ---- 1. sequence QC rows for the recovered accessions ---------------------
want <- unique(strip_version(ncbi$accession[!is.na(ncbi$accession)]))
missing_acc <- setdiff(want, unique(seq_data$accession))
log_info("{length(missing_acc)} accessions missing from sequence_qc")

if (length(missing_acc)) {
  gb_new <- fetch_gb_records(missing_acc, progress = FALSE)
  qc_new <- build_sequence_qc_table(gb_new)
  log_info("Built {nrow(qc_new)} new sequence-QC rows")
  seq_data <- dplyr::bind_rows(seq_data, qc_new) |>
    dplyr::distinct(accession, .keep_all = TRUE)
  saveRDS(seq_data, "results/sequence_qc.rds")
  log_success("sequence_qc.rds: {nrow(seq_data)} rows")
}

# ---- 2. gene annotations: backfill definitions, add recovered species -----
defs <- seq_data |>
  dplyr::filter(!is.na(definition)) |>
  dplyr::distinct(accession, definition)

if (!"definition" %in% names(genes)) genes$definition <- NA_character_

genes <- genes |>
  dplyr::mutate(acc_nov = strip_version(accession)) |>
  dplyr::rows_patch(
    defs |> dplyr::transmute(acc_nov = accession, definition),
    by = "acc_nov", unmatched = "ignore"
  )
log_info(
  "gene_annotations: {sum(!is.na(genes$definition))}/{nrow(genes)} rows now carry a definition"
)

# Species that already have at least one annotated accession need nothing;
# the rest (those recovered in step 1) get the same 5-accession subsample the
# pipeline would have drawn for them.
annotated_species <- ncbi |>
  dplyr::transmute(species = organism, acc_nov = strip_version(accession)) |>
  dplyr::filter(acc_nov %in% genes$acc_nov) |>
  dplyr::pull(species) |>
  unique()

set.seed(42)
sub_new <- ncbi |>
  dplyr::transmute(species = organism, accession) |>
  dplyr::filter(!is.na(species), !species %in% annotated_species) |>
  dplyr::group_by(species) |>
  dplyr::slice_sample(n = 5) |>
  dplyr::ungroup()
log_info("{dplyr::n_distinct(sub_new$species)} species need gene annotations")

if (nrow(sub_new)) {
  ann_new <- fetch_gene_annotations(sub_new$accession, progress = FALSE)
  ann_new$acc_nov <- strip_version(ann_new$accession)
  genes <- dplyr::bind_rows(genes, ann_new)
  log_success("Added {nrow(ann_new)} annotation rows")
}

genes <- dplyr::select(genes, -acc_nov)
saveRDS(genes, "results/genes_subsamp_50_df.rds")

# ---- 3. rebuild the derived tables ----------------------------------------
bdqc_taxo <- bdqc |>
  dplyr::select(species, vernacular_fr, vernacular_en, group_en) |>
  dplyr::distinct()
ca_risk <- load_risk_status("data/CA_especes_en_peril.csv", jurisdiction = "CA")
qc_risk <- load_risk_status("data/QC_especes_en_peril.csv", jurisdiction = "QC")

summary_table <- build_summary_dataframe(ncbi, genes, bdqc_taxo, ca_risk, qc_risk)
saveRDS(summary_table, "results/summary_table.rds")
log_success("summary_table.rds: {nrow(summary_table)} species")

species_by_accession <- ncbi |>
  dplyr::filter(!is.na(organism), !is.na(accession)) |>
  dplyr::transmute(accession = strip_version(accession), organism)

intraspecific_variation <- seq_data |>
  dplyr::filter(!is.na(seq_length), seq_length <= 10000) |>
  dplyr::mutate(marker = assign_single_marker(gene, definition)) |>
  dplyr::filter(!is.na(marker)) |>
  dplyr::inner_join(species_by_accession, by = "accession") |>
  dplyr::group_by(organism, marker) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::summarise(
    n_seq = dplyr::n(),
    median_slen = stats::median(seq_length),
    sd_slen = stats::sd(seq_length),
    cv_pct = round(sd_slen / median_slen * 100, 1),
    min_slen = min(seq_length), max_slen = max(seq_length),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(cv_pct))
saveRDS(intraspecific_variation, "results/intraspecific_variation.rds")
log_success("intraspecific_variation.rds: {nrow(intraspecific_variation)} species x marker pairs")

query_primers <- read.csv2("data/primers_map_group_bdqc_list_01122025.csv")
primers_agg <- query_primers |>
  dplyr::group_by(query_marker) |>
  dplyr::summarise(
    groups = paste(sort(unique(group)), collapse = " / "),
    markers = dplyr::first(markers), .groups = "drop"
  )
marker_stats <- ncbi |>
  dplyr::mutate(query_marker = sub(".*? AND ", "", query)) |>
  dplyr::left_join(primers_agg, by = "query_marker") |>
  dplyr::group_by(groups, markers, query_marker) |>
  dplyr::summarise(
    n_sequences = dplyr::n(), n_species = dplyr::n_distinct(organism, na.rm = TRUE),
    median_slen = stats::median(slen, na.rm = TRUE), sd_slen = stats::sd(slen, na.rm = TRUE),
    q25_slen = stats::quantile(slen, 0.25, na.rm = TRUE),
    q75_slen = stats::quantile(slen, 0.75, na.rm = TRUE), .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_sequences))
saveRDS(marker_stats, "results/marker_stats.rds")
log_success("marker_stats.rds: {nrow(marker_stats)} rows")
