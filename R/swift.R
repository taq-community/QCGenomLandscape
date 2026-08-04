#' Build the Keystone v3 password/project-scoped auth request body
#'
#' Pure helper factored out of [swift_auth()] so the request body can be
#' unit tested without performing any network call.
#'
#' @param username,password,project_name,domain Character scalars
#' @return A list, ready to be passed to [httr2::req_body_json()]
#' @noRd
build_swift_auth_body <- function(username, password, project_name, domain) {
  list(
    auth = list(
      identity = list(
        methods = list("password"),
        password = list(
          user = list(
            name = username,
            domain = list(name = domain),
            password = password
          )
        )
      ),
      scope = list(
        project = list(
          name = project_name,
          domain = list(name = domain)
        )
      )
    )
  )
}

#' Authenticate against an OpenStack Keystone endpoint
#'
#' Performs a Keystone v3 password/project-scoped auth (used by Alliance
#' Canada's Arbutus object storage, among other OpenStack Swift deployments)
#' and returns the resulting token.
#'
#' @param auth_url Character, Keystone auth URL (no trailing `/v3/...`),
#'   default `Sys.getenv("SWIFT_AUTH_URL")`
#' @param username,password,project_name Character, default from
#'   `SWIFT_USERNAME`/`SWIFT_PASSWORD`/`SWIFT_PROJECT_NAME` env vars
#' @param domain Character, Keystone domain, default `SWIFT_DOMAIN` env var
#'   or `"default"`
#' @param request_fn Function with signature `(url)` returning an `httr2`
#'   request, default [httr2::request()]
#' @return Character scalar, the Keystone auth token (from the
#'   `X-Subject-Token` response header)
#' @export
swift_auth <- function(auth_url = Sys.getenv("SWIFT_AUTH_URL"),
                        username = Sys.getenv("SWIFT_USERNAME"),
                        password = Sys.getenv("SWIFT_PASSWORD"),
                        project_name = Sys.getenv("SWIFT_PROJECT_NAME"),
                        domain = Sys.getenv("SWIFT_DOMAIN", "default"),
                        request_fn = httr2::request) {
  body <- build_swift_auth_body(username, password, project_name, domain)

  resp <- request_fn(paste0(auth_url, "/v3/auth/tokens")) |>
    httr2::req_body_json(body) |>
    httr2::req_perform()

  httr2::resp_header(resp, "X-Subject-Token")
}

#' Upload a local file to Arbutus (OpenStack Swift) object storage
#'
#' @param path Character, path to the local file to upload
#' @param object_name Character, name of the object in the container,
#'   default `basename(path)`
#' @param storage_url Character, the Swift account storage URL (not secret --
#'   this is the public endpoint shared for this project), default the
#'   QCGenomicLandscape Arbutus account
#' @param container Character, container (bucket) name, default `"QCGenomicLandscape"`
#' @param token Character, a Keystone auth token, default [swift_auth()].
#'   Pass an explicit token when uploading multiple files to avoid
#'   re-authenticating for each one.
#' @param request_fn Function with signature `(url)` returning an `httr2`
#'   request, default [httr2::request()]
#' @return Invisibly, the uploaded object's path within the container
#'   (`"{container}/{object_name}"`)
#' @export
upload_to_swift <- function(path,
                             object_name = basename(path),
                             storage_url = "https://object-arbutus.alliancecan.ca/swift/v1/cc7d995026f443268e822676773a4252",
                             container = "QCGenomicLandscape",
                             token = swift_auth(),
                             request_fn = httr2::request) {
  url <- paste(storage_url, container, object_name, sep = "/")

  request_fn(url) |>
    httr2::req_method("PUT") |>
    httr2::req_headers(`X-Auth-Token` = token) |>
    httr2::req_body_file(path) |>
    httr2::req_perform()

  invisible(paste(container, object_name, sep = "/"))
}
