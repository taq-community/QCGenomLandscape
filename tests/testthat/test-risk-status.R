test_that("load_risk_status parses CA (COSEPAC) status and drops inactive/blank rows", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines(
    c(
      "Nom.scientifique,Statut.selon.le.COSEPAC",
      "\"Alces alces (Linnaeus, 1758)\",Menacée",
      "\"Ursus americanus\",Non active",
      "\"Canis lupus\","
    ),
    tmp
  )

  result <- load_risk_status(tmp, jurisdiction = "CA")

  expect_equal(nrow(result), 1)
  expect_equal(result$species, "Alces alces")
  expect_equal(result$status, "Menacée")
  expect_equal(result$jurisdiction, "CA")
})

test_that("load_risk_status parses QC (LEMV) status and drops retired/unmonitored rows", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines(
    c(
      "GENRE,ESPECE,STATUT_LEMV",
      "Rana,pipiens,Vulnérable",
      "Ursus,americanus,Retirée",
      "Canis,lupus,Non suivie"
    ),
    tmp
  )

  result <- load_risk_status(tmp, jurisdiction = "QC")

  expect_equal(nrow(result), 1)
  expect_equal(result$species, "Rana pipiens")
  expect_equal(result$status, "Vulnérable")
})

test_that("load_risk_status translates French labels to English when translate = TRUE", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("GENRE,ESPECE,STATUT_LEMV", "Rana,pipiens,Vulnérable"), tmp)

  result <- load_risk_status(tmp, jurisdiction = "QC", translate = TRUE)

  expect_equal(result$status, "Vulnerable")
})

test_that("load_risk_status keeps French labels by default", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("GENRE,ESPECE,STATUT_LEMV", "Rana,pipiens,Vulnérable"), tmp)

  result <- load_risk_status(tmp, jurisdiction = "QC")

  expect_equal(result$status, "Vulnérable")
})
