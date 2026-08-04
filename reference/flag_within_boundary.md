# Flag points as within a polygon boundary

Flag points as within a polygon boundary

## Usage

``` r
flag_within_boundary(points_sf, boundary_sf, flag_name = "in_boundary")
```

## Arguments

- points_sf:

  `sf` points object

- boundary_sf:

  `sf` polygon object (transformed to `points_sf`'s CRS if it differs)

- flag_name:

  Character, name of the new logical column, default `"in_boundary"`

## Value

`points_sf` with an added logical column named `flag_name`
