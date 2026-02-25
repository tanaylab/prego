# Convert a MotifDB object back to a tidy data frame

Convert a MotifDB object back to a tidy data frame

## Usage

``` r
motif_db_to_dataframe(motif_db)
```

## Arguments

- motif_db:

  A MotifDB object

## Value

A tidy data frame with columns for motif, position, and nucleotide
probabilities

## Examples

``` r
head(motif_db_to_dataframe(MOTIF_DB))
#> # A tibble: 6 × 6
#>   motif                          pos       A       C       G       T
#>   <chr>                        <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 HOCOMOCO.AHR_HUMAN.H11MO.0.B     1 0.266   0.117   0.364   0.253  
#> 2 HOCOMOCO.AHR_HUMAN.H11MO.0.B     2 0.0714  0.0779  0.227   0.623  
#> 3 HOCOMOCO.AHR_HUMAN.H11MO.0.B     3 0.143   0.286   0.136   0.435  
#> 4 HOCOMOCO.AHR_HUMAN.H11MO.0.B     4 0.0195  0.00649 0.948   0.0260 
#> 5 HOCOMOCO.AHR_HUMAN.H11MO.0.B     5 0.00649 0.974   0.00649 0.0130 
#> 6 HOCOMOCO.AHR_HUMAN.H11MO.0.B     6 0.0195  0.00649 0.968   0.00649
```
