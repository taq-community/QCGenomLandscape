#' Donut chart: gene prevalence by taxonomic group
#'
#' @param genes_grouped Tibble with columns `species`, `gene_group`, `groupe`
#'   (gene classification joined to species and their finer taxonomic group,
#'   e.g. via [assign_gene_group()] and [classify_taxon_group()])
#' @param gene_groups Character vector of `gene_group` values to plot as
#'   separate facets, default `c("COI", "Photosynthesis-related (rbcL, matK, etc.)")`
#' @param min_species Integer; gene groups with fewer than this many total
#'   species (summed across taxonomic groups) are dropped, default 5
#' @return A `ggplot` object (caller is responsible for `ggplot2::ggsave()`)
#' @export
plot_gene_prevalence <- function(genes_grouped,
                                  gene_groups = c("COI", "Photosynthesis-related (rbcL, matK, etc.)"),
                                  min_species = 5) {
  gene_prevalence <- genes_grouped |>
    dplyr::distinct(species, gene_group, groupe) |>
    dplyr::group_by(groupe, gene_group) |>
    dplyr::summarise(n_sp = dplyr::n_distinct(species), .groups = "drop") |>
    dplyr::left_join(
      genes_grouped |>
        dplyr::distinct(species, groupe) |>
        dplyr::group_by(groupe) |>
        dplyr::summarise(n_total = dplyr::n_distinct(species), .groups = "drop"),
      by = "groupe"
    ) |>
    dplyr::mutate(pct = n_sp / n_total * 100)

  gene_groups_keep <- gene_prevalence |>
    dplyr::group_by(gene_group) |>
    dplyr::summarise(total = sum(n_sp), .groups = "drop") |>
    dplyr::filter(total >= min_species) |>
    dplyr::pull(gene_group)

  gene_prevalence_filtered <- gene_prevalence |>
    dplyr::filter(gene_group %in% gene_groups_keep) |>
    dplyr::filter(gene_group %in% gene_groups)

  ggplot2::ggplot(gene_prevalence_filtered, ggplot2::aes(x = stats::reorder(groupe, pct), y = pct)) +
    ggplot2::geom_col(ggplot2::aes(fill = pct), show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(pct, 1), "%")), size = 3) +
    ggplot2::coord_polar() +
    ggplot2::facet_wrap(~gene_group, ncol = 3) +
    ggplot2::scale_y_continuous(limits = c(-10, 100), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::scale_fill_gradientn(colours = rev(c("#6C5B7B", "#C06C84", "#F67280", "#F8B195"))) +
    ggplot2::labs(x = NULL, y = "Species with gene (%)", title = "Gene prevalence by taxonomic group") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 10))
}

#' Donut chart: genomic-data coverage by species-at-risk status
#'
#' @param risk_genes Tibble with columns `species`, `status`, `jurisdiction`,
#'   `has_data` (risk status joined to gene-data presence, e.g. via
#'   [load_risk_status()] joined with a `gene_group` presence flag)
#' @param jurisdiction Character, which `jurisdiction` value to filter to,
#'   default `"QC"`
#' @return A `ggplot` object (caller is responsible for `ggplot2::ggsave()`)
#' @export
plot_risk_status_coverage <- function(risk_genes, jurisdiction = "QC") {
  risk_plot <- risk_genes |>
    dplyr::filter(.data$jurisdiction == .env$jurisdiction) |>
    dplyr::group_by(status) |>
    dplyr::summarise(
      n_species = dplyr::n_distinct(species),
      n_sp_with_data = dplyr::n_distinct(species[has_data]),
      pct = n_sp_with_data / n_species * 100,
      .groups = "drop"
    )

  ggplot2::ggplot(risk_plot, ggplot2::aes(x = stats::reorder(status, pct), y = pct)) +
    ggplot2::geom_col(ggplot2::aes(fill = pct)) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(pct, 1), "%")), size = 4) +
    ggplot2::coord_polar() +
    ggplot2::scale_y_continuous(limits = c(-10, 100), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::scale_fill_gradientn("% species with data", colours = rev(c("#6C5B7B", "#C06C84", "#F67280", "#F8B195"))) +
    ggplot2::labs(
      x = NULL, y = "Species with genomic data (%)",
      title = "Genomic data coverage by species at risk status"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 14))
}
