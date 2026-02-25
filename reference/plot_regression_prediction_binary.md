# Plot the cumulative of the regression model's prediction stratified by the response variable

Plot the cumulative of the regression model's prediction stratified by
the response variable

## Usage

``` r
plot_regression_prediction_binary(pred, response)
```

## Arguments

- pred:

  the 'pred' field from the regression result

- response:

  the 'response' field from the regression result (the response
  variable). Should be binary (0/1).

## Examples

``` r
if (FALSE) { # \dontrun{
res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], score_metric = "ks")
plot_regression_prediction_binary(res_binary$pred, res_binary$response)
} # }
```
