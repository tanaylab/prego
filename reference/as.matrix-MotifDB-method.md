# Convert a MotifDB object to a matrix

Convert a MotifDB object to a matrix

## Usage

``` r
# S4 method for class 'MotifDB'
as.matrix(x, ...)
```

## Arguments

- x:

  MotifDB object

- ...:

  ignored arguments

## Value

A matrix containing the motif probabilities, rownames are
motif_position, colnames are nucleotides

## Examples

``` r
as.matrix(MOTIF_DB["HOMER.GATA3_2"])
#>                     A     C     G     T
#> HOMER.GATA3_2_1 0.662 0.066 0.006 0.266
#> HOMER.GATA3_2_2 0.001 0.007 0.991 0.001
#> HOMER.GATA3_2_3 0.989 0.004 0.001 0.006
#> HOMER.GATA3_2_4 0.002 0.023 0.001 0.974
#> HOMER.GATA3_2_5 0.825 0.061 0.011 0.103
#> HOMER.GATA3_2_6 0.778 0.048 0.129 0.045
#> HOMER.GATA3_2_7 0.184 0.401 0.348 0.067
#> HOMER.GATA3_2_8 0.433 0.167 0.359 0.041
```
