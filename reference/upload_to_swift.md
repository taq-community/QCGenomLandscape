# Upload a local file to Arbutus (OpenStack Swift) object storage

Upload a local file to Arbutus (OpenStack Swift) object storage

## Usage

``` r
upload_to_swift(
  path,
  object_name = basename(path),
  storage_url =
    "https://object-arbutus.alliancecan.ca/swift/v1/cc7d995026f443268e822676773a4252",
  container = "QCGenomicLandscape",
  token = swift_auth(),
  request_fn = httr2::request
)
```

## Arguments

- path:

  Character, path to the local file to upload

- object_name:

  Character, name of the object in the container, default
  `basename(path)`

- storage_url:

  Character, the Swift account storage URL (not secret – this is the
  public endpoint shared for this project), default the
  QCGenomicLandscape Arbutus account

- container:

  Character, container (bucket) name, default `"QCGenomicLandscape"`

- token:

  Character, a Keystone auth token, default
  [`swift_auth()`](https://taq-community.github.io/QCGenomLandscape/reference/swift_auth.md).
  Pass an explicit token when uploading multiple files to avoid
  re-authenticating for each one.

- request_fn:

  Function with signature `(url)` returning an `httr2` request, default
  [`httr2::request()`](https://httr2.r-lib.org/reference/request.html)

## Value

Invisibly, the uploaded object's path within the container
(`"{container}/{object_name}"`)
