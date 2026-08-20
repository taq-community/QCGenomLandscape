# The 16 Atlas datasets whose records are already molecular/eDNA-derived

Identified by keyword search across all ~880 Atlas dataset names
(sequence, DNA, barcode, metagenome, eDNA, INSDC, GenBank, BOLD, ITS,
etc.), then manually verified against false positives – two matches were
dropped because the keyword hit was in the *publisher's* institutional
name (e.g. "CIBIO (Research Center in Biodiversity and Genetic
Resources)"), not the dataset's actual content. Comparing eDNA against
the *whole* Atlas would be partly circular, since a meaningful share of
it is already eDNA/molecular data – some of it (iBOL) the same
underlying source as this package's own BOLD extraction.

## Usage

``` r
atlas_molecular_datasets()
```

## Value

Character vector of exact `dataset_name` values
