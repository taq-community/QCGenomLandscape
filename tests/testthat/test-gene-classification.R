test_that("assign_gene_group classifies COI variants", {
  expect_equal(assign_gene_group(c("COX1", "coi", "COXI")), rep("COI", 3))
})

test_that("assign_gene_group classifies Cytb variants", {
  expect_equal(assign_gene_group(c("cytb", "COB", "cyt b")), rep("Cytb", 3))
})

test_that("assign_gene_group classifies ND genes individually", {
  expect_equal(assign_gene_group(c("ND1", "nad1")), rep("ND1", 2))
  expect_equal(assign_gene_group(c("ND2", "nad2")), rep("ND2", 2))
  expect_equal(assign_gene_group(c("ND4", "nad4")), rep("ND4", 2))
  expect_equal(assign_gene_group(c("ND5", "nad5")), rep("ND5", 2))
})

test_that("assign_gene_group classifies rRNA genes with descriptive labels", {
  expect_equal(assign_gene_group(c("12S", "rrnS", "s-rrna")), rep("12S rRNA", 3))
  expect_equal(assign_gene_group(c("16S", "rrnL", "l-rrna")), rep("16S rRNA", 3))
  expect_equal(assign_gene_group("ITS"), "Nuclear rRNA / ITS")
})

test_that("assign_gene_group classifies photosynthesis-related genes", {
  expect_equal(assign_gene_group(c("rbcL", "matK")), rep("Photosynthesis-related (rbcL, matK, etc.)", 2))
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
