# Get the length of a MotifDB object

Get the length of a MotifDB object

## Usage

``` r
# S4 method for class 'MotifDB'
length(x)
```

## Arguments

- x:

  MotifDB object

## Value

The number of motifs in the object

## Examples

``` r
length(MOTIF_DB)
#> [1] 3867
length(MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")])
#> [1] 2
```
