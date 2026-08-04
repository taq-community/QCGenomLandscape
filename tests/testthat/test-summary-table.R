test_that("build_summary_dataframe joins sequence counts, gene markers, taxonomy, and risk status", {
  ncbi_results <- tibble::tibble(
    # a batched query (multiple species OR'd together) -- species must come
    # from `organism`, not be parsed out of `query`
    query = c(
      "(Alces alces[Organism] OR Ursus americanus[Organism]) AND COI[Gene] AND voucher[Title]",
      "(Alces alces[Organism] OR Ursus americanus[Organism]) AND COI[Gene] AND voucher[Title]"
    ),
    organism = c("Alces alces", "Alces alces"),
    accession = c("ACC1", "ACC2")
  )
  genes_df <- tibble::tibble(
    accession = c("ACC1", "ACC2"),
    gene = c("COX1", "COX1")
  )
  bdqc_taxo <- tibble::tibble(
    species = "Alces alces",
    vernacular_fr = "Orignal",
    vernacular_en = "Moose",
    group_en = "Mammals"
  )
  ca_risk <- tibble::tibble(species = "Alces alces", status = "Not at risk")
  qc_risk <- tibble::tibble(species = character(), status = character())

  result <- build_summary_dataframe(ncbi_results, genes_df, bdqc_taxo, ca_risk, qc_risk)

  expect_equal(nrow(result), 1)
  expect_equal(result$`Nom scientifique`, "Alces alces")
  expect_equal(result$`Nom commun FR`, "Orignal")
  expect_equal(result$`Groupe taxonomique`, "Mammals")
  expect_equal(result$`Séquences totales (NCBI)`, 2)
  expect_equal(result$COI, 2)
  expect_equal(result$`Statut de l'espèce au Canada`, "Not at risk")
  expect_true(is.na(result$`Statut de l'espèce au Québec`))
})

test_that("build_summary_dataframe drops rows with no species (empty input)", {
  ncbi_results <- tibble::tibble(query = character(0), organism = character(0), accession = character(0))
  genes_df <- tibble::tibble(accession = character(0), gene = character(0))
  bdqc_taxo <- tibble::tibble(
    species = character(0), vernacular_fr = character(0),
    vernacular_en = character(0), group_en = character(0)
  )
  ca_risk <- tibble::tibble(species = character(0), status = character(0))
  qc_risk <- tibble::tibble(species = character(0), status = character(0))

  result <- build_summary_dataframe(ncbi_results, genes_df, bdqc_taxo, ca_risk, qc_risk)

  expect_equal(nrow(result), 0)
})
