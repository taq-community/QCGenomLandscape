#' Classify BDQC taxonomic groups into finer display categories
#'
#' Reclassifies the coarse `group_en` column of the BDQC species list into
#' finer taxonomic categories used for figures (e.g. splitting mammals into
#' marine/terrestrial, arthropods into insects/crustaceans).
#'
#' @param group_en Character vector, BDQC's coarse taxonomic group
#' @param order Character vector, taxonomic order (same length as `group_en`)
#' @param class Character vector, taxonomic class (same length as `group_en`)
#' @param phylum Character vector, taxonomic phylum (same length as `group_en`)
#' @param kingdom Character vector, taxonomic kingdom (same length as `group_en`)
#' @return Character vector of finer taxonomic group labels
#' @examples
#' classify_taxon_group("Mammals", order = "Cetacea")
#' classify_taxon_group("Arthropods", class = "Insecta")
#' @export
classify_taxon_group <- function(group_en, order = NA, class = NA, phylum = NA, kingdom = NA) {
  dplyr::case_when(
    group_en == "Fish" ~ "Fish",
    group_en == "Mammals" & order %in% c("Cetacea", "Pinnipedia", "Sirenia") ~ "Marine mammals",
    group_en == "Mammals" ~ "Terrestrial mammals",
    group_en == "Amphibians" ~ "Amphibians",
    group_en == "Reptiles" ~ "Reptiles",
    group_en == "Birds" ~ "Birds",
    group_en == "Arthropods" & class == "Insecta" ~ "Insects",
    group_en == "Arthropods" &
      class %in% c("Malacostraca", "Branchiopoda", "Copepoda", "Maxillopoda", "Ostracoda") ~ "Crustaceans",
    group_en == "Other invertebrates" & phylum == "Mollusca" ~ "Mollusks",
    group_en %in% c("Angiosperms", "Conifers", "Bryophytes",
                     "Vascular cryptogam", "Other plants", "Algae") ~ "Plants",
    group_en == "Fungi" ~ "Fungi",
    group_en == "Other taxons" & kingdom == "Bacteria" ~ "Bacteria",
    group_en == "Other taxons" & kingdom == "Protozoa" ~ "Protozoa",
    TRUE ~ "Other"
  )
}
