#' Classify a gene name into a coarse marker group
#'
#' Unifies the two gene-classification rule sets that previously diverged
#' between `create_dataframe.R` (short labels, `Cytb`/`ND1`/`ND2`/`ND4`/`ND5`)
#' and `taxon_representation.R` (descriptive `rRNA`/photosynthesis labels).
#'
#' A `;`-joined list of several gene names (e.g. all 13 mitochondrial
#' protein-coding genes annotated together on one GenBank record) is
#' classified as `"Multi-gene / genome-scale record"` before any of the
#' single-marker patterns are tried -- a full mitogenome isn't comparable to
#' a short single-marker fragment (can't be meaningfully aligned against
#' one, wildly different expected length, etc.), so grouping it with actual
#' single-COI-marker records under a shared label would be misleading for
#' anything downstream that treats `gene_group` as "same kind of sequence"
#' (e.g. [flag_length_outliers()], [flag_barcode_gap_outliers()]).
#'
#' @param gene Character vector of raw gene names (e.g. `"COX1"`, `"cytb"`)
#' @return Character vector, one of `"COI"`, `"COII"`, `"COIII"`, `"Cytb"`,
#'   `"ND1"`, `"ND2"`, `"ND3"`, `"ND4"`, `"ND4L"`, `"ND5"`, `"ND6"`,
#'   `"ATP6"`, `"ATP8"`, `"12S rRNA"`, `"16S rRNA"`, `"Nuclear rRNA / ITS"`,
#'   `"Fungal protein-coding (RPB1/RPB2, TEF1)"`,
#'   `"Photosynthesis-related (rbcL, matK, etc.)"`,
#'   `"Multi-gene / genome-scale record"`, `"Other"`, or `NA` if `gene` is `NA`
#' @examples
#' assign_gene_group(c("COX1", "cytb", "12S", "rbcL", "trnL", "xyz123"))
#' assign_gene_group("ND1;ND2;COX1;COX2;ATP8;ATP6;COX3;ND3;ND4L;ND4;ND5;ND6;CYTB")
#' @export
assign_gene_group <- function(gene) {
  gene_lower <- tolower(gene)
  dplyr::case_when(
    is.na(gene_lower) ~ NA_character_,
    grepl(";", gene_lower, fixed = TRUE) ~ "Multi-gene / genome-scale record",
    grepl("^(cox1|coi|coxi)$", gene_lower) ~ "COI",
    grepl("^(cox2|coii)$", gene_lower) ~ "COII",
    grepl("^(cox3|coiii)$", gene_lower) ~ "COIII",
    grepl("^(cytb|cob|cyt b|cytochrome b)$", gene_lower) ~ "Cytb",
    grepl("^(nd1|nad1)$", gene_lower) ~ "ND1",
    grepl("^(nd2|nad2)$", gene_lower) ~ "ND2",
    grepl("^(nd3|nad3)$", gene_lower) ~ "ND3",
    grepl("^(nd4l|nad4l)$", gene_lower) ~ "ND4L",
    grepl("^(nd4|nad4)$", gene_lower) ~ "ND4",
    grepl("^(nd5|nad5)$", gene_lower) ~ "ND5",
    grepl("^(nd6|nad6)$", gene_lower) ~ "ND6",
    grepl("^atp6$", gene_lower) ~ "ATP6",
    grepl("^atp8$", gene_lower) ~ "ATP8",
    grepl("(12s|rrns|s-rrna)", gene_lower) ~ "12S rRNA",
    grepl("(16s|rrnl|l-rrna)", gene_lower) ~ "16S rRNA",
    grepl("^(18s|28s|25s|5\\.8s|5s|its[12]?|its|ssu|lsu)( rrna)?$", gene_lower) ~ "Nuclear rRNA / ITS",
    grepl("^(rpb1|rpb2|tef1)$", gene_lower) ~ "Fungal protein-coding (RPB1/RPB2, TEF1)",
    grepl("^(rbcl|matk|psba|ndhf|trnh-psba|trnl)", gene_lower) ~ "Photosynthesis-related (rbcL, matK, etc.)",
    TRUE ~ "Other"
  )
}
