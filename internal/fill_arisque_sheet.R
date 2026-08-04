library(dplyr)
library(stringr)
library(googlesheets4)
library(taxadb)

SHEET_ID <- "1W5__6vlrnY1IeZPYEeTAx7a8bG8oLTpj4GRISthXsos"

# Severity order: index 1 = most severe
COSEPAC_ORDER <- c(
  "Disparue du pays",
  "En voie de disparition",
  "Menacée",
  "Préoccupante",
  "Données insuffisantes"
)

LEMV_ORDER <- c(
  "Disparue",
  "En voie de disparition",
  "Menacée",
  "Vulnérable",
  "Susceptible d'être désignée",
  "Candidate"
)

highest_status <- function(statuses, order) {
  matches <- order[order %in% statuses]
  if (length(matches) == 0) return(NA_character_)
  matches[1]
}

# -----------------------------------------------------------------------------
# BDQC taxonomic reference
# -----------------------------------------------------------------------------

bdqc <- read.csv("data/bdqc_list_01122025.csv") |>
  filter(rank == "species") |>
  select(species, class, order, family, phylum, group_en) |>
  distinct(species, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# qPCR kits
# -----------------------------------------------------------------------------

kits <- read.csv("data/qPCR_kits.csv", check.names = FALSE) |>
  mutate(species = str_extract(`Scientific name`, "^\\w+\\s+\\w+")) |>
  group_by(species) |>
  summarise(kit_ref = paste(sort(unique(Reference)), collapse = "; "), .groups = "drop")

# -----------------------------------------------------------------------------
# Canada (LEP / COSEPAC) — une ligne par espèce, statut le plus sévère
# -----------------------------------------------------------------------------

ca_data <- read.csv("data/CA_especes_en_peril.csv", fileEncoding = "UTF-8-BOM") |>
  mutate(
    species        = str_extract(`Nom.scientifique`, "^\\w+\\s+\\w+"),
    nom_commun_fr  = `Nom.commun.utilisé.par.le.COSEPAC`,
    statut_cosepac = `Statut.selon.le.COSEPAC`,
    statut_lep     = `Statut.à.l.annexe.1`,
    type           = `Groupe.taxinomique`
  ) |>
  filter(!statut_cosepac %in% c("", "Non active", "Non en péril")) |>
  group_by(species) |>
  summarise(
    nom_commun_fr  = first(nom_commun_fr),
    type           = first(type),
    statut_cosepac = highest_status(statut_cosepac, COSEPAC_ORDER),
    statut_lep     = highest_status(statut_lep[statut_lep != ""], COSEPAC_ORDER),
    .groups        = "drop"
  )

# -----------------------------------------------------------------------------
# Québec (LEMV) — une ligne par espèce, statut le plus sévère
# -----------------------------------------------------------------------------

qc_data <- read.csv("data/QC_especes_en_peril.csv", fileEncoding = "UTF-8-BOM") |>
  mutate(species = paste(GENRE, ESPECE)) |>
  filter(!STATUT_LEMV %in% c("Retirée", "Non suivie", "Aucun statut")) |>
  group_by(species, CLASSE, ORDRE, FAMILLE, Nom_francais, Nom_anglais, GRAND_GROUPE) |>
  summarise(
    statut_lemv = highest_status(STATUT_LEMV, LEMV_ORDER),
    .groups     = "drop"
  ) |>
  rename(
    nom_commun_fr = Nom_francais,
    nom_commun_en = Nom_anglais,
    classe        = CLASSE,
    ordre         = ORDRE,
    famille       = FAMILLE,
    type_qc       = GRAND_GROUPE
  )

# -----------------------------------------------------------------------------
# taxize fallback — fill missing class/order/family/phylum via NCBI
# -----------------------------------------------------------------------------

base_df <- full_join(ca_data, qc_data, by = "species") |>
  left_join(bdqc, by = "species") |>
  mutate(
    classe  = coalesce(classe, class),
    ordre   = coalesce(ordre,  order),
    famille = coalesce(famille, family),
    phylum  = phylum
  )

missing_taxo <- base_df |>
  filter(is.na(classe) | is.na(ordre) | is.na(famille)) |>
  pull(species) |>
  unique()

if (length(missing_taxo) > 0) {
  message("Looking up ", length(missing_taxo), " species via taxadb (ITIS)...")

  # td_create("itis") must be run once to download the local database
  taxadb_df <- taxadb::filter_name(missing_taxo, provider = "itis") |>
    select(
      species    = scientificName,
      classe_tz  = class,
      ordre_tz   = order,
      famille_tz = family,
      phylum_tz  = phylum
    ) |>
    distinct(species, .keep_all = TRUE)

  base_df <- base_df |>
    left_join(taxadb_df, by = "species") |>
    mutate(
      classe  = coalesce(classe,  classe_tz),
      ordre   = coalesce(ordre,   ordre_tz),
      famille = coalesce(famille, famille_tz),
      phylum  = coalesce(phylum,  phylum_tz)
    ) |>
    select(-classe_tz, -ordre_tz, -famille_tz, -phylum_tz)
}

# -----------------------------------------------------------------------------
# NCBI full genome results
# -----------------------------------------------------------------------------

genomes <- tryCatch(
  readRDS("results/ncbi_full_genome_species_at_risk.rds") |>
    mutate(
      nuclear_urls = ifelse(nuclear_genome,  "Oui", "Non"),
      mito_urls    = ifelse(mitochondrial_genome, "Oui", "Non")
    ) |>
    select(species, nuclear_urls, mito_urls),
  error = function(e) {
    message("ncbi_full_genome_species_at_risk.rds not found — genome columns will be empty.")
    tibble(species = character(), nuclear_urls = character(), mito_urls = character())
  }
)

# -----------------------------------------------------------------------------
# Full join — une ligne par espèce, outer pour les sans-match
# -----------------------------------------------------------------------------

arisque <- base_df |>
  mutate(
    `Nom commun` = coalesce(nom_commun_fr.x, nom_commun_fr.y),
    `Type` = case_when(
      # Aquatic invertebrates — crustaceans, molluscs, other aquatic invertebrates
      classe %in% c("Malacostraca", "Branchiopoda", "Copepoda", "Maxillopoda",
                    "Ostracoda", "Bivalvia", "Gastropoda", "Cephalopoda",
                    "Polychaeta", "Hirudinea", "Oligochaeta") ~ "Invertébrés aquatiques",
      phylum %in% c("Mollusca", "Annelida", "Echinodermata",
                    "Cnidaria", "Porifera")                   ~ "Invertébrés aquatiques",
      # Aquatic insects (orders with fully aquatic larval stages)
      classe == "Insecta" & ordre %in% c("Ephemeroptera", "Plecoptera",
                                          "Trichoptera", "Odonata",
                                          "Megaloptera")      ~ "Invertébrés aquatiques",
      # Terrestrial invertebrates
      classe %in% c("Insecta", "Arachnida", "Chilopoda",
                    "Diplopoda", "Collembola")                ~ "Invertébrés terrestres",
      group_en %in% c("Arthropods", "Other invertebrates")   ~ "Invertébrés terrestres",
      # Keep original type for all other groups
      TRUE ~ coalesce(type, type_qc)
    )
  ) |>
  left_join(kits, by = "species") |>
  left_join(genomes, by = "species") |>
  transmute(
    `Type`                          = `Type`,
    `Classe`                        = classe,
    `Ordre`                         = ordre,
    `Famille`                       = famille,
    `Espèce`                        = species,
    `Nom commun`                    = `Nom commun`,
    `Common Name`                   = nom_commun_en,
    `Statut COSEPAC`                = statut_cosepac,
    `Statut LEP (Annexe 1)`         = statut_lep,
    `Statut LEMV`                   = statut_lemv,
    `Kit disponible`                = kit_ref,
    `Full nuclear genome sequenced?`    = nuclear_urls,
    `Full mitogenome sequenced?`        = mito_urls
  ) |>
  arrange(`Type`, `Classe`, `Ordre`, `Famille`, `Espèce`)


# -----------------------------------------------------------------------------
# Write to Google Sheets — "À risque" tab
# -----------------------------------------------------------------------------

gs4_auth()

sheet_write(
  data  = arisque,
  ss    = SHEET_ID,
  sheet = "À risque"
)

message("Done — ", nrow(arisque), " rows written to 'À risque'.")
