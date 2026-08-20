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

Character vector, one of `"COI"`, `"COII"`, `"COIII"`, `"Cytb"`,
`"ND1"`, `"ND2"`, `"ND3"`, `"ND4"`, `"ND4L"`, `"ND5"`, `"ND6"`,
`"ATP6"`, `"ATP8"`, `"12S rRNA"`, `"16S rRNA"`, `"Nuclear rRNA / ITS"`,
`"Fungal protein-coding (RPB1/RPB2, TEF1)"`,
`"Photosynthesis-related (rbcL, matK, etc.)"`,
`"Multi-gene / genome-scale record"`, `"Other"`, or `NA` if `gene` is
`NA`

## Details

A `;`-joined list of several gene names (e.g. all 13 mitochondrial
protein-coding genes annotated together on one GenBank record) is
classified as `"Multi-gene / genome-scale record"` before any of the
single-marker patterns are tried – a full mitogenome isn't comparable to
a short single-marker fragment (can't be meaningfully aligned against
one, wildly different expected length, etc.), so grouping it with actual
single-COI-marker records under a shared label would be misleading for
anything downstream that treats `gene_group` as "same kind of sequence"
(e.g.
[`flag_length_outliers()`](https://taq-community.github.io/QCGenomLandscape/reference/flag_length_outliers.md),
[`flag_barcode_gap_outliers()`](https://taq-community.github.io/QCGenomLandscape/reference/flag_barcode_gap_outliers.md)).

## Examples

``` r
assign_gene_group(c("COX1", "cytb", "12S", "rbcL", "trnL", "xyz123"))
#> [1] "COI"                                      
#> [2] "Cytb"                                     
#> [3] "12S rRNA"                                 
#> [4] "Photosynthesis-related (rbcL, matK, etc.)"
#> [5] "Photosynthesis-related (rbcL, matK, etc.)"
#> [6] "Other"                                    
assign_gene_group("ND1;ND2;COX1;COX2;ATP8;ATP6;COX3;ND3;ND4L;ND4;ND5;ND6;CYTB")
#> [1] "Multi-gene / genome-scale record"
```
