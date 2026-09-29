test_that("classify_taxon_group splits mammals into marine/terrestrial", {
  expect_equal(classify_taxon_group("Mammals", order = "Cetacea"), "Marine mammals")
  expect_equal(classify_taxon_group("Mammals", order = "Pinnipedia"), "Marine mammals")
  expect_equal(classify_taxon_group("Mammals", order = "Carnivora"), "Terrestrial mammals")
})

test_that("classify_taxon_group splits arthropods into insects/crustaceans", {
  expect_equal(classify_taxon_group("Arthropods", class = "Insecta"), "Insects")
  expect_equal(classify_taxon_group("Arthropods", class = "Malacostraca"), "Crustaceans")
})

test_that("classify_taxon_group handles mollusks via phylum", {
  expect_equal(classify_taxon_group("Other invertebrates", phylum = "Mollusca"), "Mollusks")
})

test_that("classify_taxon_group handles bacteria/protozoa via kingdom", {
  expect_equal(classify_taxon_group("Other taxons", kingdom = "Bacteria"), "Bacteria")
  expect_equal(classify_taxon_group("Other taxons", kingdom = "Protozoa"), "Protozoa")
})

test_that("classify_taxon_group falls back to Other for unmapped groups", {
  expect_equal(classify_taxon_group("Unmapped Group"), "Other")
})

test_that("classify_taxon_group is vectorized", {
  result <- classify_taxon_group(
    group_en = c("Fish", "Birds", "Fungi"),
    order = c(NA, NA, NA)
  )
  expect_equal(result, c("Fish", "Birds", "Fungi"))
})

test_that("classify_taxon_group gives arachnids their own group", {
  expect_equal(classify_taxon_group("Arthropods", class = "Arachnida"), "Arachnids")
  # still distinct from the insect/crustacean branches above it
  expect_equal(classify_taxon_group("Arthropods", class = "Collembola"), "Other")
})

test_that("classify_taxon_group routes diatoms and brown algae to Algae", {
  # filed under "Other taxons" in the backbone, so group_en alone misses them
  expect_equal(
    classify_taxon_group("Other taxons", class = "Bacillariophyceae",
                         phylum = "Ochrophyta", kingdom = "Chromista"),
    "Algae"
  )
  expect_equal(
    classify_taxon_group("Other taxons", class = "Phaeophyceae",
                         phylum = "Ochrophyta", kingdom = "Chromista"),
    "Algae"
  )
  expect_equal(classify_taxon_group("Algae"), "Algae")
})

test_that("classify_taxon_group keeps water moulds out of Algae", {
  # Oomycota are fungus-like chromists, not algae
  expect_equal(
    classify_taxon_group("Other taxons", phylum = "Oomycota", kingdom = "Chromista"),
    "Other"
  )
})

test_that("classify_taxon_group no longer folds algae into Plants", {
  expect_equal(classify_taxon_group("Angiosperms"), "Plants")
  expect_equal(classify_taxon_group("Bryophytes"), "Plants")
  expect_false(classify_taxon_group("Algae") == "Plants")
})
