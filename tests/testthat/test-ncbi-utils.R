test_that("build_organism_clause OR-combines species names, wrapped in parens", {
  expect_equal(
    build_organism_clause(c("Alces alces", "Ursus americanus")),
    "(Alces alces[Organism] OR Ursus americanus[Organism])"
  )
})

test_that("build_organism_clause wraps even a single species in parens", {
  expect_equal(build_organism_clause("Alces alces"), "(Alces alces[Organism])")
})

test_that("chunk_species splits a vector into batch_size-sized chunks, in order", {
  result <- chunk_species(c("A", "B", "C", "D", "E"), batch_size = 2)

  expect_length(result, 3)
  expect_equal(result[[1]], c("A", "B"))
  expect_equal(result[[2]], c("C", "D"))
  expect_equal(result[[3]], "E")
})

test_that("chunk_species returns one chunk when batch_size exceeds the vector length", {
  result <- chunk_species(c("A", "B"), batch_size = 25)

  expect_length(result, 1)
  expect_equal(result[[1]], c("A", "B"))
})
