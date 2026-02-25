# Export a motif regression model

Export a motif regression model

## Usage

``` r
export_regression_model(model, fn = NULL)
```

## Arguments

- model:

  a motif regression model, as returned by `regress_pwm` with
  `motif_num = 1`

- fn:

  a file name to save the model to. If NULL - the model is returned as a
  list

## Value

None

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7
)
export_fn <- tempfile()
export_regression_model(export_fn)
r <- load_regression(export_fn)
} # }
```
