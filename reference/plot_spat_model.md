# Plot spatial model of the regression result

Plot spatial model of the regression result

## Usage

``` r
plot_spat_model(spat, title = "Spatial model")
```

## Arguments

- spat:

  the 'spat' field from the regression result

- title:

  a title for the plot (optional)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(sequences_example, response_mat_example)
plot_spat_model(res$spat)
} # }
```
