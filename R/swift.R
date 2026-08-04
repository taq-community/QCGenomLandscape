#' Build the Keystone v3 application-credential auth request body
#'
#' Pure helper factored out of [swift_auth()] so the request body can be
#' unit tested without performing any network call.
#'
#' Application-credential auth needs no username/password/project/domain --
#' the credential ID already uniquely identifies the user and is pre-scoped
#' to a project when created. This is also what makes it work with MFA-
#' protected accounts: the credential is generated once (interactively,
#' after an MFA challenge) via the Horizon dashboard or CLI, and from then on
#' authenticates on its own.
#'
#' @param app_cred_id,app_cred_secret Character scalars
#' @return A list, ready to be passed to [httr2::req_body_json()]
#' @noRd
build_swift_auth_body <- function(app_cred_id, app_cred_secret) {
  list(
    auth = list(
      identity = list(
        methods = list("application_credential"),
        application_credential = list(
          id = app_cred_id,
          secret = app_cred_secret
        )
      )
    )
  )
}

#' Authenticate against an OpenStack Keystone endpoint
#'
#' Performs a Keystone v3 application-credential auth (used by Alliance
#' Canada's Arbutus object storage, among other OpenStack Swift deployments)
#' and returns the resulting token. Application credentials -- not username/
#' password -- because Arbutus accounts protected by MFA can't authenticate
#' with a plain password over the API; an application credential is
#' generated once interactively (Horizon dashboard: Identity > Application
#' Credentials, or `openstack application credential create <name>` after an
#' interactive/MFA login) and used from then on.
#'
#' @param auth_url Character, Keystone auth URL (trailing slash tolerated),
#'   default `SWIFT_AUTH_URL` env var or Arbutus's identity endpoint (not
#'   secret -- specific to this project, but fine to override for another one)
#' @param app_cred_id,app_cred_secret Character, default `SWIFT_APP_CRED_ID`/
#'   `SWIFT_APP_CRED_SECRET` env vars (no default -- these are secret and
#'   must be set)
#' @param request_fn Function with signature `(url)` returning an `httr2`
#'   request, default [httr2::request()]
#' @return Character scalar, the Keystone auth token (from the
#'   `X-Subject-Token` response header)
#' @export
swift_auth <- function(auth_url = Sys.getenv("SWIFT_AUTH_URL", "https://identity.arbutus.alliancecan.ca"),
                        app_cred_id = Sys.getenv("SWIFT_APP_CRED_ID"),
                        app_cred_secret = Sys.getenv("SWIFT_APP_CRED_SECRET"),
                        request_fn = httr2::request) {
  auth_url <- sub("/+$", "", auth_url)
  body <- build_swift_auth_body(app_cred_id, app_cred_secret)

  logger::log_info("Authenticating against Keystone ({auth_url})...")

  resp <- request_fn(paste0(auth_url, "/v3/auth/tokens")) |>
    httr2::req_body_json(body) |>
    httr2::req_perform()

  logger::log_success("Keystone auth OK")

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
  object_size <- file.info(path)$size

  logger::log_info("Uploading {path} ({format(object_size, big.mark = ',')} bytes) -> {container}/{object_name}")

  tryCatch(
    {
      request_fn(url) |>
        httr2::req_method("PUT") |>
        httr2::req_headers(`X-Auth-Token` = token) |>
        httr2::req_body_file(path) |>
        httr2::req_perform()
    },
    error = function(e) {
      logger::log_error("Swift upload failed for {object_name}: {e$message}")
      stop(e)
    }
  )

  logger::log_success("Uploaded {container}/{object_name}")

  invisible(paste(container, object_name, sep = "/"))
}
