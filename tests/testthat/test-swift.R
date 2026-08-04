test_that("build_swift_auth_body builds a Keystone v3 application-credential body", {
  body <- build_swift_auth_body(
    app_cred_id = "abc123",
    app_cred_secret = "s3cret"
  )

  expect_equal(body$auth$identity$methods, list("application_credential"))
  expect_equal(body$auth$identity$application_credential$id, "abc123")
  expect_equal(body$auth$identity$application_credential$secret, "s3cret")
})

test_that("build_swift_auth_body needs no username/password/project/domain fields", {
  # Application-credential auth is pre-scoped to a project at creation time --
  # there should be no separate `scope`/`password` block in the request body.
  body <- build_swift_auth_body("abc123", "s3cret")

  expect_null(body$auth$scope)
  expect_null(body$auth$identity$password)
})

test_that("build_swift_auth_body never leaks the secret outside the secret field", {
  body <- build_swift_auth_body("abc123", "s3cret")

  expect_false(grepl("s3cret", body$auth$identity$application_credential$id, fixed = TRUE))
})
