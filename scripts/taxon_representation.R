library(ggplot2)
bb_taxo <- read.csv("data/bdqc_list_01122025.csv") |>
    dplyr::filter(rank == "species")

ncbi_results <- readRDS("results/ncbi_results.rds")
ncbi_results$species <- stringr::str_extract(ncbi_results$query, "^\\w+\\s+\\w+")

length(unique(ncbi_results$species)) / length(unique(bb_taxo$species))

count_by_phylum <- ncbi_results |>
    dplyr::left_join(
        bb_taxo |> dplyr::select(species, phylum) |> dplyr::distinct(),
        by = "species"
    ) |>
    dplyr::group_by(phylum) |>
    dplyr::summarise(count = dplyr::n()) |>
    dplyr::mutate(source = "ncbi") |>
    dplyr::bind_rows(
        bb_taxo |>
            dplyr::group_by(phylum) |>
            dplyr::summarise(count = sum(obs_count)) |>
            dplyr::mutate(source = "biodiversitéQuébec")
    ) |>
    dplyr::filter(!is.na(phylum))

count_by_phylum |>
    dplyr::group_by(phylum) |>
    dplyr::filter(sum(count) > 500) |>
    dplyr::ungroup() |>
    dplyr::group_by(source) |>
    dplyr::mutate(proportion = count / sum(count) * 100) |>
    dplyr::ungroup() |>
    ggplot2::ggplot(ggplot2::aes(y = reorder(phylum, proportion), x = proportion, fill = source)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
        y = "Phylum",
        x = "Percentage (%)",
        title = "Relative Representation by Phylum across Sources",
        fill = "Source"
    ) +
    ggplot2::theme_light(base_size = 26) 




# Contribution spatial


# Contribution temporelle


# Contribution taxonomique


# Par type de milieu
# 1. Lac
# 2. Forêt

1. Metabar sans assignation, partir plus de la localisation
2. Pour une séquence de référence, Cb de séquence on a qui correspond à cette séquence


