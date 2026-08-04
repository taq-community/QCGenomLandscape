# Authenticate against an OpenStack Keystone endpoint

Performs a Keystone v3 application-credential auth (used by Alliance
Canada's Arbutus object storage, among other OpenStack Swift
deployments) and returns the resulting token. Application credentials –
not username/ password – because Arbutus accounts protected by MFA can't
authenticate with a plain password over the API; an application
credential is generated once interactively (Horizon dashboard: Identity
\> Application Credentials, or
`openstack application credential create <name>` after an
interactive/MFA login) and used from then on.

## Usage

``` r
swift_auth(
  auth_url = Sys.getenv("SWIFT_AUTH_URL", "https://identity.arbutus.alliancecan.ca"),
  app_cred_id = Sys.getenv("SWIFT_APP_CRED_ID"),
  app_cred_secret = Sys.getenv("SWIFT_APP_CRED_SECRET"),
  request_fn = httr2::request
)
```

## Arguments

- auth_url:

  Character, Keystone auth URL (trailing slash tolerated), default
  `SWIFT_AUTH_URL` env var or Arbutus's identity endpoint (not secret –
  specific to this project, but fine to override for another one)

- app_cred_id, app_cred_secret:

  Character, default `SWIFT_APP_CRED_ID`/ `SWIFT_APP_CRED_SECRET` env
  vars (no default – these are secret and must be set)

- request_fn:

  Function with signature `(url)` returning an `httr2` request, default
  [`httr2::request()`](https://httr2.r-lib.org/reference/request.html)

## Value

Character scalar, the Keystone auth token (from the `X-Subject-Token`
response header)
