# Plot response variable averages vs the regression model's prediction

this would return a scatter plot of the response variable averages vs
the regression model's prediction

## Usage

``` r
plot_regression_prediction(pred, response, point_size = 0.5, alpha = 1)
```

## Arguments

- pred:

  the 'pred' field from the regression result

- response:

  the 'response' field from the regression result (the response
  variable)

- point_size:

  the size of the points in the plot (default: 0.5)

- alpha:

  the transparency of the points in the plot (default: 1)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(sequences_example, response_mat_example)
plot_regression_prediction(res$pred, res$response)
} # }
```
