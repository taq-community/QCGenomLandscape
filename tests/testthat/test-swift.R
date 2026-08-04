test_that("build_swift_auth_body builds a Keystone v3 password/project-scoped body", {
  body <- build_swift_auth_body(
    username = "svissault",
    password = "s3cret",
    project_name = "def-someproj",
    domain = "default"
  )

  expect_equal(body$auth$identity$methods, list("password"))
  expect_equal(body$auth$identity$password$user$name, "svissault")
  expect_equal(body$auth$identity$password$user$password, "s3cret")
  expect_equal(body$auth$identity$password$user$domain$name, "default")
  expect_equal(body$auth$scope$project$name, "def-someproj")
  expect_equal(body$auth$scope$project$domain$name, "default")
})

test_that("build_swift_auth_body never leaks the password outside the password field", {
  body <- build_swift_auth_body("user", "s3cret", "proj", "default")

  expect_false(grepl("s3cret", body$auth$scope$project$name, fixed = TRUE))
})
