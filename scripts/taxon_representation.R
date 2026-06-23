library(ggplot2)
bb_taxo <- read.csv("data/bdqc_list_01122025.csv") |>
    dplyr::filter(rank == "species")

# Reclassify group_en into finer taxonomic categories
bb_taxo$groupe <- dplyr::case_when(
    bb_taxo$group_en == "Fish" ~ "Fish",
    bb_taxo$group_en == "Mammals" &
        bb_taxo$order %in% c("Cetacea", "Pinnipedia", "Sirenia") ~ "Marine mammals",
    bb_taxo$group_en == "Mammals" ~ "Terrestrial mammals",
    bb_taxo$group_en == "Amphibians" ~ "Amphibians",
    bb_taxo$group_en == "Reptiles" ~ "Reptiles",
    bb_taxo$group_en == "Birds" ~ "Birds",
    bb_taxo$group_en == "Arthropods" &
        bb_taxo$class == "Insecta" ~ "Insects",
    bb_taxo$group_en == "Arthropods" &
        bb_taxo$class %in% c("Malacostraca", "Branchiopoda", "Copepoda",
                             "Maxillopoda", "Ostracoda") ~ "Crustaceans",
    bb_taxo$group_en == "Other invertebrates" &
        bb_taxo$phylum == "Mollusca" ~ "Mollusks",
    bb_taxo$group_en %in% c("Angiosperms", "Conifers", "Bryophytes",
                            "Vascular cryptogam", "Other plants", "Algae") ~ "Plants",
    bb_taxo$group_en == "Fungi" ~ "Fungi",
    bb_taxo$group_en == "Other taxons" &
        bb_taxo$kingdom == "Bacteria" ~ "Bacteria",
    bb_taxo$group_en == "Other taxons" &
        bb_taxo$kingdom == "Protozoa" ~ "Protozoa",
    TRUE ~ "Other"
)

ncbi_results <- readRDS("results/ncbi_results.rds")
ncbi_results$species <- stringr::str_extract(ncbi_results$query, "^\\w+\\s+\\w+")

# Accession number by species
library(ggplot2)
library(dplyr)

set.seed(42)
accessions <- ncbi_results |>
    group_by(species) |>
    slice_sample(n = 5) |>
    ungroup() |>
    pull(accession)

rentrez::set_entrez_key(Sys.getenv("NCBI_API_KEY"))
batches <- split(accessions, ceiling(seq_along(accessions) / 200))
genes_df <- purrr::map_df(seq_along(batches), function(i) {
    batch <- batches[[i]]
    tryCatch({
        gb_raw <- rentrez::entrez_fetch(db = "nucleotide", id = batch, rettype = "gb", retmode = "xml")
        gb <- xml2::read_xml(gb_raw, options = "HUGE")
        seq_nodes <- xml2::xml_find_all(gb, ".//GBSeq")
        purrr::map_df(seq_nodes, function(seq_node) {
            acc <- xml2::xml_text(xml2::xml_find_first(seq_node, ".//GBSeq_accession-version"))
            gene_nodes <- xml2::xml_find_all(seq_node, ".//GBFeature[GBFeature_key='gene']")
            purrr::map_df(gene_nodes, function(node) {
                gene_name <- xml2::xml_text(
                    xml2::xml_find_first(node, ".//GBQualifier[GBQualifier_name='gene']/GBQualifier_value")
                )
                location <- xml2::xml_text(xml2::xml_find_first(node, ".//GBFeature_location"))
                dplyr::tibble(accession = acc, gene = gene_name, location = location)
            })
        })
    }, error = function(e) {
        warning(sprintf("Batch %d/%d failed: %s", i, length(batches), e$message))
        dplyr::tibble(accession = character(), gene = character(), location = character())
    })
}, .progress = TRUE)

##
saveRDS(genes_df, "results/genes_subsamp_50_df.rds")
library(dplyr)
genes_df <- readRDS("results/genes_subsamp_50_df.rds")


## Overall
data_fig <- ncbi_results |>
    filter(accession %in% genes_df$accession) |>
    nest_join(genes_df, by = "accession") |>
    mutate(has_gene = is.null(genes_df)) |>
    full_join(bb_taxo, by = "species")

plot_df <- data_fig |>
    mutate(has_gene = sapply(genes_df, \(x) !is.null(x) && nrow(x) > 0)) |>
    filter(!is.na(groupe)) |>
    group_by(groupe) |>
    summarise(
        n_species = n_distinct(species),
        n_sp_with_genes = n_distinct(species[has_gene]),
        pct = n_sp_with_genes / n_species * 100
    )

ggplot(plot_df, aes(x = reorder(groupe, pct), y = pct)) +
    geom_col(aes(fill = pct)) +
    geom_text(aes(label = paste0(round(pct, 1), "%")), size = 6) +
    coord_polar() +
    labs(
        x = NULL,
        y = "Species with genomic data (%)",
        title = "Proportion of species with gene data by group"
    ) +
    # Scale y axis so bars don't start in the center
    scale_y_continuous(
        limits = c(-10, 100),
        expand = c(0, 0),
        breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_fill_gradientn(
        "% species with data",
        colours = rev(c("#6C5B7B", "#C06C84", "#F67280", "#F8B195"))
    ) +
    theme_minimal(base_size = 24)


# Gene grouping function
assign_gene_group <- function(gene) {
    gene_lower <- tolower(gene)
    dplyr::case_when(
        # COI / COX1
        grepl("^(cox1|coi|coxi)$", gene_lower) ~ "COI",
        # Cytochrome b
        # grepl("^(cytb|cob|cyt b|cytochrome b)$", gene_lower) ~ "Cytochrome b",
        # Other COX subunits
        # grepl("^(cox2|coii|cox3|coiii|coxii|coxiii)$", gene_lower) ~ "COX2/COX3",
        # NADH dehydrogenase
        # grepl("^(nd[1-6]|nd4l|nad[1-6]|nad4l)$", gene_lower) ~ "NADH dehydrogenase (ND)",
        # ATP synthase
        # grepl("^(atp[68]|atpase[68])$", gene_lower) ~ "ATP synthase (ATP6/8)",
        # Ribosomal RNA (mitochondrial)
        grepl("(12s|rrns|s-rrna)", gene_lower) ~ "12S rRNA",
        grepl("(16s|rrnl|l-rrna)", gene_lower) ~ "16S rRNA",
        grepl("^(18s|28s|5\\.8s|5s|its[12]?|its)$", gene_lower) ~ "Nuclear rRNA / ITS",
        # tRNA
        # grepl("^trn", gene_lower) ~ "tRNA",
        # Ribosomal proteins (plastid/mito)
        # grepl("^(rps|rpl)", gene_lower) ~ "Ribosomal proteins",
        # Photosynthesis-related (rbcL, matK, etc.)
        grepl("^(rbcl|matk|psba|ndhf|trnh-psba)", gene_lower) ~ "Photosynthesis-related (rbcL, matK, etc.)",
        # Control region / D-loop
        # grepl("(d-loop|control region|cr)", gene_lower) ~ "Control region",
        TRUE ~ "Other"
    )
}

# Build gene-level prevalence per taxonomic group
genes_grouped <- genes_df |>
    mutate(gene_group = assign_gene_group(gene)) |>
    left_join(
        ncbi_results |> select(accession, species),
        by = "accession"
    ) |>
    full_join(
        bb_taxo |> select(species, groupe) |> distinct(),
        by = "species"
    ) |>
    filter(!is.na(groupe))

gene_prevalence <- genes_grouped |>
    distinct(species, gene_group, groupe) |>
    group_by(groupe, gene_group) |>
    summarise(n_sp = n_distinct(species), .groups = "drop") |>
    left_join(
        genes_grouped |>
            distinct(species, groupe) |>
            group_by(groupe) |>
            summarise(n_total = n_distinct(species), .groups = "drop"),
        by = "groupe"
    ) |>
    mutate(pct = n_sp / n_total * 100)

# Keep only gene groups with enough data to be meaningful
gene_groups_keep <- gene_prevalence |>
    group_by(gene_group) |>
    summarise(total = sum(n_sp)) |>
    filter(total >= 5) |>
    pull(gene_group)

gene_prevalence_filtered <- gene_prevalence |>
    filter(gene_group %in% gene_groups_keep) |>
    filter(gene_group %in% c("COI", "Photosynthesis-related (rbcL, matK, etc.)"))
    # filter(!is.na(gene_group) | gene_group != "Other")

ggplot(gene_prevalence_filtered, aes(x = reorder(groupe, pct), y = pct)) +
    geom_col(aes(fill = pct), show.legend = FALSE) +
    geom_text(aes(label = paste0(round(pct, 1), "%")), size = 3) +
    coord_polar() +
    facet_wrap(~ gene_group, ncol = 3) +
    scale_y_continuous(
        limits = c(-10, 100),
        expand = c(0, 0),
        breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_fill_gradientn(
        colours = rev(c("#6C5B7B", "#C06C84", "#F67280", "#F8B195"))
    ) +
    labs(
        x = NULL,
        y = "Species with gene (%)",
        title = "Gene prevalence by taxonomic group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        strip.text = element_text(face = "bold", size = 10)
    )

ggsave("results/genes_prevalence.svg", width = 18, height = 18)

## Species at risk — genomic coverage by status
ca_risk <- read.csv("data/CA_especes_en_peril.csv", fileEncoding = "UTF-8-BOM") |>
    mutate(
        species = stringr::str_extract(`Nom.scientifique`, "^\\w+\\s+\\w+"),
        status = `Statut.selon.le.COSEPAC`,
        jurisdiction = "CA"
    ) |>
    filter(status != "" & status != "Non active") |>
    select(species, status, jurisdiction) |>
    distinct()

qc_risk <- read.csv("data/QC_especes_en_peril.csv", fileEncoding = "UTF-8-BOM") |>
    mutate(
        species = paste(GENRE, ESPECE),
        status = STATUT_LEMV,
        jurisdiction = "QC"
    ) |>
    filter(status != "Retirée" & status != "Non suivie") |>
    select(species, status, jurisdiction) |>
    distinct()

risk_df <- bind_rows(ca_risk, qc_risk) |>
    mutate(status = dplyr::recode(status,
        "Disparue" = "Extinct",
        "Disparue du pays" = "Extirpated",
        "Données insuffisantes" = "Data deficient",
        "En voie de disparition" = "Endangered",
        "Menacée" = "Threatened",
        "Non en péril" = "Not at risk",
        "Préoccupante" = "Special concern",
        "Candidate" = "Candidate",
        "Susceptible" = "Likely to be designated",
        "Vulnérable" = "Vulnerable"
    ))

# Join with gene data
risk_genes <- risk_df |>
    left_join(
        genes_grouped |>
            distinct(species, gene_group),
        by = "species"
    ) |>
    mutate(has_data = !is.na(gene_group))

risk_plot <- risk_genes |>
    filter(jurisdiction == "QC") |>
    group_by(status) |>
    summarise(
        n_species = n_distinct(species),
        n_sp_with_data = n_distinct(species[has_data]),
        pct = n_sp_with_data / n_species * 100,
        .groups = "drop"
    )

ggplot(risk_plot, aes(x = reorder(status, pct), y = pct)) +
    geom_col(aes(fill = pct)) +
    geom_text(aes(label = paste0(round(pct, 1), "%")), size = 4) +
    coord_polar() +
    scale_y_continuous(
        limits = c(-10, 100),
        expand = c(0, 0),
        breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_fill_gradientn(
        "% species with data",
        colours = rev(c("#6C5B7B", "#C06C84", "#F67280", "#F8B195"))
    ) +
    labs(
        x = NULL,
        y = "Species with genomic data (%)",
        title = "Genomic data coverage by species at risk status"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        strip.text = element_text(face = "bold", size = 14)
    )

ggsave("results/risk_status_coverage.svg", width = 16, height = 10)
