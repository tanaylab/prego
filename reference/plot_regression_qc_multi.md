# Plot the regression results for multiple motifs

plot the regression results when `motif_num` \> 1

## Usage

``` r
plot_regression_qc_multi(
  reg,
  title = glue("Motif regression results (consensus: {reg$consensus})"),
  subtitle = NULL,
  caption = NULL,
  point_size = 0.01,
  alpha = 0.5,
  response = NULL
)
```

## Arguments

- reg:

  output of `regress_pwm`

- title:

  a title for the plot (optional)

- subtitle:

  a subtitle for the plot (optional)

- caption:

  a caption for the plot (optional). When caption is NULL a default
  caption would be plotted.

- point_size:

  the size of the points in the scatter plot

- alpha:

  the transparency of the points in the scatter plot

- response:

  the response variable

## Examples

``` r
if (FALSE) { # \dontrun{
res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 3], motif_num = 3)
plot_regression_qc_multi(res_binary)
} # }
```
