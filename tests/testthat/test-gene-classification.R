test_that("assign_gene_group classifies COI variants", {
  expect_equal(assign_gene_group(c("COX1", "coi", "COXI")), rep("COI", 3))
})

test_that("assign_gene_group classifies Cytb variants", {
  expect_equal(assign_gene_group(c("cytb", "COB", "cyt b")), rep("Cytb", 3))
})

test_that("assign_gene_group classifies COII/COIII variants", {
  expect_equal(assign_gene_group(c("COX2", "coii")), rep("COII", 2))
  expect_equal(assign_gene_group(c("COX3", "coiii")), rep("COIII", 2))
})

test_that("assign_gene_group classifies ND genes individually", {
  expect_equal(assign_gene_group(c("ND1", "nad1")), rep("ND1", 2))
  expect_equal(assign_gene_group(c("ND2", "nad2")), rep("ND2", 2))
  expect_equal(assign_gene_group(c("ND3", "nad3")), rep("ND3", 2))
  expect_equal(assign_gene_group(c("ND4L", "nad4l")), rep("ND4L", 2))
  expect_equal(assign_gene_group(c("ND4", "nad4")), rep("ND4", 2))
  expect_equal(assign_gene_group(c("ND5", "nad5")), rep("ND5", 2))
  expect_equal(assign_gene_group(c("ND6", "nad6")), rep("ND6", 2))
})

test_that("assign_gene_group classifies ATP synthase subunit genes", {
  expect_equal(assign_gene_group("ATP6"), "ATP6")
  expect_equal(assign_gene_group("ATP8"), "ATP8")
})

test_that("assign_gene_group classifies rRNA genes with descriptive labels", {
  expect_equal(assign_gene_group(c("12S", "rrnS", "s-rrna")), rep("12S rRNA", 3))
  expect_equal(assign_gene_group(c("16S", "rrnL", "l-rrna")), rep("16S rRNA", 3))
  expect_equal(assign_gene_group("ITS"), "Nuclear rRNA / ITS")
  # observed in real GenBank gene fields (not just the bare "5.8S"/"SSU" forms)
  expect_equal(assign_gene_group(c("5.8S rRNA", "25S rRNA", "SSU", "LSU")), rep("Nuclear rRNA / ITS", 4))
})

test_that("assign_gene_group classifies fungal protein-coding markers", {
  expect_equal(assign_gene_group(c("RPB1", "RPB2", "TEF1")), rep("Fungal protein-coding (RPB1/RPB2, TEF1)", 3))
})

test_that("assign_gene_group classifies photosynthesis-related genes", {
  expect_equal(assign_gene_group(c("rbcL", "matK")), rep("Photosynthesis-related (rbcL, matK, etc.)", 2))
  # trnL: chloroplast marker in the same barcoding family as rbcL/matK
  expect_equal(assign_gene_group("trnL"), "Photosynthesis-related (rbcL, matK, etc.)")
})

test_that("assign_gene_group classifies a ';'-joined multi-gene record separately from any single marker", {
  multi <- "ND1;ND2;COX1;COX2;ATP8;ATP6;COX3;ND3;ND4L;ND4;ND5;ND6;CYTB"
  expect_equal(assign_gene_group(multi), "Multi-gene / genome-scale record")
  # even when every individual gene listed would otherwise resolve cleanly
  expect_equal(assign_gene_group("COX1;COX1"), "Multi-gene / genome-scale record")
})

test_that("assign_gene_group is case-insensitive", {
  expect_equal(assign_gene_group("Cox1"), "COI")
})

test_that("assign_gene_group falls back to Other for unmatched genes", {
  expect_equal(assign_gene_group("xyz123"), "Other")
})

test_that("assign_gene_group returns NA for NA input", {
  expect_equal(assign_gene_group(NA_character_), NA_character_)
})

test_that("assign_single_marker splits the composite gene groups", {
  # the whole point: these three share one assign_gene_group() label but are
  # not length-comparable, so a per-marker metric must keep them apart
  expect_equal(
    assign_single_marker(c("rbcL", "matK", "trnL")),
    c("RBCL", "MATK", "TRNL")
  )
  expect_equal(assign_single_marker(c("RPB1", "RPB2")), c("RPB1", "RPB2"))
  # case variants and the optional " rRNA" suffix are the same marker
  expect_equal(assign_single_marker(c("matK", "matk", "MATK")), rep("MATK", 3))
  expect_equal(assign_single_marker(c("18S rRNA", "18S")), c("18S", "18S"))
})

test_that("assign_single_marker keeps synonym unification where the label is one marker", {
  expect_equal(assign_single_marker(c("COX1", "COI")), c("COI", "COI"))
  expect_equal(assign_single_marker(c("cytb", "cob")), c("Cytb", "Cytb"))
})

test_that("assign_single_marker returns NA where no single marker applies", {
  expect_true(is.na(assign_single_marker("ND1;ND2;COX1"))) # genome-scale
  expect_true(is.na(assign_single_marker("xyz123")))       # "Other"
  expect_true(is.na(assign_single_marker(NA_character_)))
})

test_that("assign_gene_group falls back to the DEFINITION line when /gene= is absent", {
  # the fungal barcode: ITS regions are annotated as misc_RNA/rRNA with
  # /product=, so `gene` is NA and these records used to vanish entirely
  defs <- c(
    "Abortiporus biennis strain CBS 148.40 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence",
    "Acarospora fuscata isolate OTU1995 5.8S ribosomal RNA gene, partial sequence; internal transcribed spacer 2, complete sequence"
  )
  expect_equal(
    assign_gene_group(rep(NA_character_, 2), defs),
    rep("Nuclear rRNA / ITS", 2)
  )
})

test_that("the DEFINITION fallback never overrides an annotated gene", {
  expect_equal(
    assign_gene_group("COX1", "Homo sapiens internal transcribed spacer 1"),
    "COI"
  )
})

test_that("assign_gene_group is unchanged when no definition is supplied", {
  expect_true(is.na(assign_gene_group(NA_character_)))
})

test_that("classify_definition_group recognises the other marker families", {
  expect_equal(
    classify_definition_group("Lepomis gibbosus cytochrome b (cytb) gene, partial cds"),
    "Cytb"
  )
  expect_equal(
    classify_definition_group("Danaus plexippus cytochrome c oxidase subunit I gene, partial cds"),
    "COI"
  )
  expect_equal(
    classify_definition_group("Boreogadus saida mitochondrion, complete genome"),
    "Multi-gene / genome-scale record"
  )
})

test_that("assign_single_marker labels definition-derived ITS records", {
  expect_equal(
    assign_single_marker(NA_character_, "strain X internal transcribed spacer 1, partial sequence"),
    "ITS"
  )
})
