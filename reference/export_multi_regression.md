# Export a multiple motif regression model

Export a multiple motif regression model

## Usage

``` r
export_multi_regression(reg, fn = NULL)
```

## Arguments

- reg:

  a multiple motif regression model, as returned by `regress_pwm` with
  `motif_num > 1`

- fn:

  a file name to save the model to. If NULL - the model is returned as a
  list

## Value

None

## Examples

``` r
if (FALSE) { # \dontrun{
res_multi <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7,
    motif_num = 2
)
export_fn <- tempfile()
export_multi_regression(res_multi, export_fn)

light_res <- export_multi_regression(res_multi)

# loading can be done by:
r <- load_multi_regression(export_fn)
} # }
```
