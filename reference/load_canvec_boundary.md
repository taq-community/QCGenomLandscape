# Load a CanVec administrative-boundary layer, optionally filtered

Load a CanVec administrative-boundary layer, optionally filtered

## Usage

``` r
load_canvec_boundary(
  gdb_path,
  layer = "geo_political_region_2",
  jurisdiction = NULL,
  country = NULL,
  crs = 4326
)
```

## Arguments

- gdb_path:

  Character, path to the `.gdb` directory (e.g.
  `"data/canvec_1M_CA_Admin.gdb"`)

- layer:

  Character, layer name, default `"geo_political_region_2"`

- jurisdiction:

  Integer or `NULL`; if set, keeps rows where
  `jurisdiction == jurisdiction` (e.g. `102` for Quebec)

- country:

  Integer or `NULL`; if set, keeps rows where `country == country` (e.g.
  `140` for Canada)

- crs:

  Target CRS, default `4326`

## Value

An `sf` object
