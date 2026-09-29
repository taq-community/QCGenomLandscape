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
#' Records carrying no `/gene=` tag at all fall back to their GenBank
#' DEFINITION line when `definition` is supplied. This matters more than it
#' sounds: the ITS region -- the fungal barcode -- is annotated as `misc_RNA`
#' and `rRNA` features with `/product=` qualifiers, almost never with
#' `/gene=`, so 147,740 of the 150,849 untagged records in the 2026-08-19
#' corpus are ITS. Classifying on `gene` alone discards every one of them and
#' reports fungi as ~4% covered when the true figure is an order of magnitude
#' higher.
#'
#' @param gene Character vector of raw gene names (e.g. `"COX1"`, `"cytb"`)
#' @param definition Character vector of GenBank DEFINITION lines, same
#'   length as `gene` (e.g. [parse_gb_records()]'s `definition` column).
#'   Used only where `gene` is `NA`. Default `NULL` keeps the historical
#'   gene-only behaviour.
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
assign_gene_group <- function(gene, definition = NULL) {
  gene_lower <- tolower(gene)
  group <- dplyr::case_when(
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

  if (is.null(definition)) {
    return(group)
  }
  if (length(definition) != length(gene)) {
    stop("`definition` must be the same length as `gene`", call. = FALSE)
  }
  untagged <- is.na(group) & !is.na(definition)
  group[untagged] <- classify_definition_group(definition[untagged])
  group
}

#' Classify a GenBank DEFINITION line into a coarse marker group
#'
#' Fallback for records with no `/gene=` qualifier, used by
#' [assign_gene_group()]. Order matters: a single ITS deposit routinely names
#' the small subunit, ITS1, 5.8S, ITS2 and the large subunit in one
#' DEFINITION, so the nuclear-rRNA test has to come before the narrower
#' ribosomal tests or those records scatter across several groups.
#'
#' @param definition Character vector of GenBank DEFINITION lines
#' @return Character vector of [assign_gene_group()] labels
#' @examples
#' classify_definition_group(
#'   "Laccaria amethystina internal transcribed spacer 1, partial sequence"
#' )
#' @export
classify_definition_group <- function(definition) {
  d <- tolower(definition)
  dplyr::case_when(
    is.na(d) ~ NA_character_,
    grepl("complete genome|complete sequence, whole genome|chromosome", d) ~ "Multi-gene / genome-scale record",
    grepl("internal transcribed spacer|\\bits[12]?\\b|5\\.8s|18s ribosomal|28s ribosomal|small subunit ribosomal|large subunit ribosomal", d) ~ "Nuclear rRNA / ITS",
    grepl("cytochrome c? ?oxidase subunit (i|1)\\b", d) ~ "COI",
    grepl("cytochrome b\\b", d) ~ "Cytb",
    grepl("12s ribosomal", d) ~ "12S rRNA",
    grepl("16s ribosomal", d) ~ "16S rRNA",
    grepl("\\brbcl\\b|ribulose-1,5|maturase k|\\bmatk\\b|\\btrnl\\b", d) ~ "Photosynthesis-related (rbcL, matK, etc.)",
    grepl("\\brpb[12]\\b|translation elongation factor", d) ~ "Fungal protein-coding (RPB1/RPB2, TEF1)",
    TRUE ~ "Other"
  )
}

#' Resolve a gene annotation to a single comparable marker
#'
#' [assign_gene_group()] is deliberately coarse: three of its labels are
#' *buckets* of several distinct markers -- `"Photosynthesis-related"` pools
#' rbcL (~1400 bp), matK (~850 bp) and trnL (~300 bp), `"Nuclear rRNA / ITS"`
#' pools 5.8S (~160 bp) with LSU (~3.5 kb), and `"Fungal protein-coding"`
#' pools RPB1 with RPB2. That is fine for coverage figures ("does this species
#' have a plastid marker at all?") but wrong for any length-based quality
#' metric, which needs sequences that are actually comparable to each other.
#'
#' This resolves those buckets down to the raw gene name (normalised for case
#' and for the optional ` rRNA` suffix, so `"18S rRNA"` and `"18S"` are one
#' marker), and keeps [assign_gene_group()]'s label everywhere else -- that is
#' where the label already denotes exactly one marker *and* usefully unifies
#' synonyms (`COX1`/`COI`, `cytb`/`cob`).
#'
#' Genome-scale records and the `"Other"` catch-all resolve to `NA`: neither
#' is a single marker, so neither belongs in a per-marker comparison.
#'
#' @param gene Character vector of raw gene names, as annotated on the record
#' @param definition Character vector of GenBank DEFINITION lines, same
#'   length as `gene`; passed through to [assign_gene_group()] so that
#'   records with no `/gene=` tag (overwhelmingly ITS) resolve to `"ITS"`
#'   instead of dropping out of every per-marker statistic. Default `NULL`.
#' @return Character vector of marker labels, `NA` where no single marker applies
#' @examples
#' assign_single_marker(c("COX1", "rbcL", "matK", "18S rRNA", "18S"))
#' assign_single_marker("ND1;ND2;COX1") # NA -- genome-scale record
#' @export
assign_single_marker <- function(gene, definition = NULL) {
  group <- assign_gene_group(gene, definition)
  composite <- c(
    "Photosynthesis-related (rbcL, matK, etc.)",
    "Nuclear rRNA / ITS",
    "Fungal protein-coding (RPB1/RPB2, TEF1)"
  )
  raw <- toupper(sub("\\s*rRNA$", "", trimws(gene), ignore.case = TRUE))
  # Records classified from their DEFINITION have no gene name to resolve
  # down to; an ITS deposit spanning ITS1+5.8S+ITS2 is one marker in
  # practice and is labelled as such rather than dropped for want of a tag.
  if (!is.null(definition)) {
    untagged <- is.na(gene) & !is.na(group)
    raw[untagged] <- ifelse(
      grepl("internal transcribed spacer|\\bits[12]?\\b", tolower(definition[untagged])),
      "ITS", NA_character_
    )
  }
  dplyr::case_when(
    is.na(group) ~ NA_character_,
    group %in% c("Multi-gene / genome-scale record", "Other") ~ NA_character_,
    group %in% composite ~ raw,
    TRUE ~ group
  )
}
