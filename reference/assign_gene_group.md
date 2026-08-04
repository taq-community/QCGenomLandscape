# Classify a gene name into a coarse marker group

Unifies the two gene-classification rule sets that previously diverged
between `create_dataframe.R` (short labels,
`Cytb`/`ND1`/`ND2`/`ND4`/`ND5`) and `taxon_representation.R`
(descriptive `rRNA`/photosynthesis labels).

## Usage

``` r
assign_gene_group(gene)
```

## Arguments

- gene:

  Character vector of raw gene names (e.g. `"COX1"`, `"cytb"`)

## Value

Character vector, one of `"COI"`, `"Cytb"`, `"ND1"`, `"ND2"`, `"ND4"`,
`"ND5"`, `"12S rRNA"`, `"16S rRNA"`, `"Nuclear rRNA / ITS"`,
`"Photosynthesis-related (rbcL, matK, etc.)"`, `"Other"`, or `NA` if
`gene` is `NA`

## Examples

``` r
assign_gene_group(c("COX1", "cytb", "12S", "rbcL", "xyz123"))
#> [1] "COI"                                      
#> [2] "Cytb"                                     
#> [3] "12S rRNA"                                 
#> [4] "Photosynthesis-related (rbcL, matK, etc.)"
#> [5] "Other"                                    
```
