% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-regression.R
\name{plot_regression_prediction_binary}
\alias{plot_regression_prediction_binary}
\title{Plot the cumulative of the regression model's prediction stratified by the response variable}
\usage{
plot_regression_prediction_binary(pred, response)
}
\arguments{
\item{pred}{the 'pred' field from the regression result}

\item{response}{the 'response' field from the regression result (the response variable). Should be binary (0/1).}
}
\description{
Plot the cumulative of the regression model's prediction stratified by the response variable
}
\examples{
\dontrun{
res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], score_metric = "ks")
plot_regression_prediction_binary(res_binary$pred, res_binary$response)
}

}
