# Load a motif regression model from a file

Load a motif regression model from a file

## Usage

``` r
load_regression_model(fn)
```

## Arguments

- fn:

  file name or a list with the model

## Value

a list with the following elements:

- pssm::

  a data frame with the PSSM

- spat::

  a data frame with the spatial profile

- spat_min::

  a numeric value with the minimum value of the spatial profile

- spat_max::

  a numeric value with the maximum value of the spatial profile

- bidirect::

  a boolean value indicating whether the model is bidirectional

- seq_length::

  a numeric value with the length of the sequences

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
