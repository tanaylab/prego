# Convert a MotifDB object to a data frame

Convert a MotifDB object to a data frame

## Usage

``` r
# S4 method for class 'MotifDB'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A MotifDB object

- row.names:

  NULL or a character vector giving the row names for the data frame

- optional:

  logical. If TRUE, setting row names and converting column names (to
  syntactic names: see make.names) is optional

- ...:

  additional arguments to be passed to or from methods

## Value

A data frame containing the motif probabilities

## Examples

``` r
dataset <- as.data.frame(MOTIF_DB)
head(dataset)
#>                          motif pos           A           C           G
#> 1 HOCOMOCO.AHR_HUMAN.H11MO.0.B   1 0.266233766 0.116883117 0.363636364
#> 2 HOCOMOCO.AHR_HUMAN.H11MO.0.B   2 0.071428571 0.077922078 0.227272727
#> 3 HOCOMOCO.AHR_HUMAN.H11MO.0.B   3 0.142857143 0.285714286 0.136363636
#> 4 HOCOMOCO.AHR_HUMAN.H11MO.0.B   4 0.019480519 0.006493506 0.948051948
#> 5 HOCOMOCO.AHR_HUMAN.H11MO.0.B   5 0.006493506 0.974025974 0.006493506
#> 6 HOCOMOCO.AHR_HUMAN.H11MO.0.B   6 0.019480519 0.006493506 0.967532468
#>             T
#> 1 0.253246753
#> 2 0.623376623
#> 3 0.435064935
#> 4 0.025974026
#> 5 0.012987013
#> 6 0.006493506
nrow(dataset)
#> [1] 49269
length(unique(dataset$motif))
#> [1] 3867
```
