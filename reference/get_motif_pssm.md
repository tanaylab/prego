# Extract pssm of sequences from a motif database

Extract pssm of sequences from a motif database

## Usage

``` r
get_motif_pssm(motif, dataset = all_motif_datasets())
```

## Arguments

- motif:

  name of the motif to extract from the dataset

- dataset:

  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
  additional column 'motif' containing the motif name, for example
  `HOMER_motifs` or `JASPAR_motifs`, or
  [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md).

## Value

a data frame with the pssm of the motif

## Examples

``` r
get_motif_pssm("JASPAR.HNF1A")
#> # A tibble: 15 × 5
#>      pos       A        C        G        T
#>    <dbl>   <dbl>    <dbl>    <dbl>    <dbl>
#>  1     0 0.368   0.173    0.247    0.213   
#>  2     1 0.459   0.0273   0.472    0.0418  
#>  3     2 0.0228  0.0643   0.0548   0.858   
#>  4     3 0.147   0.110    0.0275   0.715   
#>  5     4 0.982   0.000208 0.0171   0.000971
#>  6     5 0.934   0.0418   0.000528 0.0240  
#>  7     6 0.0804  0.0267   0.000442 0.893   
#>  8     7 0.236   0.265    0.321    0.178   
#>  9     8 0.934   0.000594 0.0224   0.0430  
#> 10     9 0.0452  0.000827 0.0538   0.900   
#> 11    10 0.00291 0.0158   0.000277 0.981   
#> 12    11 0.815   0.0137   0.0770   0.0943  
#> 13    12 0.878   0.0431   0.0680   0.0108  
#> 14    13 0.0865  0.822    0.0381   0.0531  
#> 15    14 0.281   0.242    0.161    0.316   
```
