#' Classify a gene name into a coarse marker group
#'
#' Unifies the two gene-classification rule sets that previously diverged
#' between `create_dataframe.R` (short labels, `Cytb`/`ND1`/`ND2`/`ND4`/`ND5`)
#' and `taxon_representation.R` (descriptive `rRNA`/photosynthesis labels).
#'
#' @param gene Character vector of raw gene names (e.g. `"COX1"`, `"cytb"`)
#' @return Character vector, one of `"COI"`, `"Cytb"`, `"ND1"`, `"ND2"`,
#'   `"ND4"`, `"ND5"`, `"12S rRNA"`, `"16S rRNA"`, `"Nuclear rRNA / ITS"`,
#'   `"Photosynthesis-related (rbcL, matK, etc.)"`, `"Other"`, or `NA` if
#'   `gene` is `NA`
#' @examples
#' assign_gene_group(c("COX1", "cytb", "12S", "rbcL", "xyz123"))
#' @export
assign_gene_group <- function(gene) {
  gene_lower <- tolower(gene)
  dplyr::case_when(
    is.na(gene_lower) ~ NA_character_,
    grepl("^(cox1|coi|coxi)$", gene_lower) ~ "COI",
    grepl("^(cytb|cob|cyt b|cytochrome b)$", gene_lower) ~ "Cytb",
    grepl("^(nd1|nad1)$", gene_lower) ~ "ND1",
    grepl("^(nd2|nad2)$", gene_lower) ~ "ND2",
    grepl("^(nd4|nad4)$", gene_lower) ~ "ND4",
    grepl("^(nd5|nad5)$", gene_lower) ~ "ND5",
    grepl("(12s|rrns|s-rrna)", gene_lower) ~ "12S rRNA",
    grepl("(16s|rrnl|l-rrna)", gene_lower) ~ "16S rRNA",
    grepl("^(18s|28s|5\\.8s|5s|its[12]?|its)$", gene_lower) ~ "Nuclear rRNA / ITS",
    grepl("^(rbcl|matk|psba|ndhf|trnh-psba)", gene_lower) ~ "Photosynthesis-related (rbcL, matK, etc.)",
    TRUE ~ "Other"
  )
}
