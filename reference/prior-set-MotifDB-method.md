# Set a new prior for a MotifDB object

Set a new prior for a MotifDB object

## Usage

``` r
# S4 method for class 'MotifDB'
prior(object) <- value
```

## Arguments

- object:

  MotifDB object

- value:

  New prior value between 0 and 1

## Value

Updated MotifDB object with new prior

## Examples

``` r
prior(MOTIF_DB)
#> [1] 0.01
prior(MOTIF_DB) <- 0.2
prior(MOTIF_DB)
#> [1] 0.2
```
