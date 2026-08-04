# Build a hexagonal summary grid of point counts within a boundary

Build a hexagonal summary grid of point counts within a boundary

## Usage

``` r
make_hex_grid(
  points_sf,
  boundary_sf,
  cellsize = 10000,
  crs_proj = 3347,
  log_scale = TRUE
)
```

## Arguments

- points_sf:

  `sf` points object (any CRS)

- boundary_sf:

  `sf` polygon object defining the grid extent

- cellsize:

  Numeric, grid cell size in meters (in `crs_proj`), default 10000

- crs_proj:

  Numeric EPSG code for the projection used to build the grid, default
  `3347` (Canada Lambert Conformal Conic)

- log_scale:

  Logical; if `TRUE`, adds a `log_count = log(n_records + 1)` column

## Value

`sf` polygon grid (WGS84), filtered to cells with `n_records > 0`
