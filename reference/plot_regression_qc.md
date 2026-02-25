# Plot the regression results

Plot QC of the regression results

## Usage

``` r
plot_regression_qc(
  reg,
  response = NULL,
  title = glue("Motif regression results (consensus: {reg$consensus})"),
  subtitle = NULL,
  caption = NULL,
  point_size = 0.5,
  alpha = 0.5
)
```

## Arguments

- reg:

  output of `regress_pwm`

- response:

  the response variable

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

## Value

a patchwork object

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(sequences_example, response_mat_example)
plot_regression_qc(res)

res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], screen_db = TRUE)
plot_regression_qc(res_binary)
} # }
```
