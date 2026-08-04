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
