library(dplyr)
library(stringr)

# Load source data
ncbi_results <- readRDS("results/ncbi_results.rds")
ncbi_results$species <- str_extract(ncbi_results$query, "^\\w+\\s+\\w+")

genes_df <- readRDS("results/genes_subsamp_50_df.rds")

bb_taxo <- read.csv("data/bdqc_list_01122025.csv") |>
    filter(rank == "species") |>
    select(species, vernacular_fr, vernacular_en, group_en) |>
    distinct()

ca_risk <- read.csv("data/CA_especes_en_peril.csv", fileEncoding = "UTF-8-BOM") |>
    mutate(
        species = str_extract(`Nom.scientifique`, "^\\w+\\s+\\w+"),
        statut_canada = `Statut.selon.le.COSEPAC`
    ) |>
    filter(statut_canada != "" & statut_canada != "Non active") |>
    select(species, statut_canada) |>
    distinct()

qc_risk <- read.csv("data/QC_especes_en_peril.csv", fileEncoding = "UTF-8-BOM") |>
    mutate(
        species = paste(GENRE, ESPECE),
        statut_quebec = STATUT_LEMV
    ) |>
    filter(statut_quebec != "Retirée" & statut_quebec != "Non suivie") |>
    select(species, statut_quebec) |>
    distinct()

# Classify genes (reusing logic from taxon_representation.R)
assign_gene_group <- function(gene) {
    gene_lower <- tolower(gene)
    dplyr::case_when(
        grepl("^(cox1|coi|coxi)$", gene_lower) ~ "COI",
        grepl("^(cytb|cob|cyt b|cytochrome b)$", gene_lower) ~ "Cytb",
        grepl("^(nd1|nad1)$", gene_lower) ~ "ND1",
        grepl("^(nd2|nad2)$", gene_lower) ~ "ND2",
        grepl("^(nd4|nad4)$", gene_lower) ~ "ND4",
        grepl("^(nd5|nad5)$", gene_lower) ~ "ND5",
        grepl("(12s|rrns|s-rrna)", gene_lower) ~ "12S",
        grepl("(16s|rrnl|l-rrna)", gene_lower) ~ "16S",
        TRUE ~ "Other"
    )
}

# Count sequences per species per gene
gene_counts <- genes_df |>
  mutate(gene_group = assign_gene_group(gene)) |>
  left_join(ncbi_results |> select(accession, species), by = "accession") |>
  filter(!is.na(species), gene_group != "Other") |>
  distinct(species, accession, gene_group) |>
  count(species, gene_group, name = "n_seq") |>
  tidyr::pivot_wider(names_from = gene_group, values_from = n_seq, values_fill = 0) |>
  distinct()

# Total sequences per species (all genes combined)
total_seq <- ncbi_results |>
    count(species, name = "n_total_seq")

# Build final summary data.frame
summary_df <- ncbi_results |>
    distinct(species) |>
    filter(!is.na(species)) |>
    left_join(bb_taxo, by = "species") |>
    left_join(ca_risk, by = "species") |>
    left_join(qc_risk, by = "species") |>
    left_join(total_seq, by = "species") |>
    left_join(gene_counts, by = "species") |>
    select(
        `Nom scientifique` = species,
        `Nom commun FR` = vernacular_fr,
        `Nom commun EN` = vernacular_en,
        `Groupe taxonomique` = group_en,
        `Statut de l'espèce au Canada` = statut_canada,
        `Statut de l'espèce au Québec` = statut_quebec,
        `Séquences totales (NCBI)` = n_total_seq,
        any_of(c("COI", "Cytb", "ND1", "ND2", "ND4", "ND5", "12S", "16S"))
    )

summary_df
