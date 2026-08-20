# Build a full hexagonal grid over a boundary, with stable cell IDs

Unlike
[`make_hex_grid()`](https://taq-community.github.io/QCGenomLandscape/reference/make_hex_grid.md)
(which counts one point layer and drops empty cells), this returns
*every* cell with a stable `cell_id` regardless of occupancy, so
multiple point layers (e.g. eDNA vs traditional occurrences) can be
joined against the same reusable grid and compared cell-for-cell.

## Usage

``` r
build_hex_grid_cells(boundary_sf, cellsize = 10000, crs_proj = 3347)
```

## Arguments

- boundary_sf:

  `sf` polygon object defining the grid extent

- cellsize:

  Numeric, grid cell size in meters (in `crs_proj`), default 10000

- crs_proj:

  Numeric EPSG code for the projection used to build the grid, default
  `3347` (Canada Lambert Conformal Conic)

## Value

`sf` polygon grid (WGS84) with a `cell_id` column
