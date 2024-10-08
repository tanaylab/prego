% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-regression.R
\name{plot_regression_prediction}
\alias{plot_regression_prediction}
\title{Plot response variable averages vs the regression model's prediction}
\usage{
plot_regression_prediction(pred, response, point_size = 0.5, alpha = 1)
}
\arguments{
\item{pred}{the 'pred' field from the regression result}

\item{response}{the 'response' field from the regression result (the response variable)}

\item{point_size}{the size of the points in the plot (default: 0.5)}

\item{alpha}{the transparency of the points in the plot (default: 1)}
}
\description{
this would return a scatter plot of the response variable averages vs the regression model's prediction
}
\examples{
\dontrun{
res <- regress_pwm(sequences_example, response_mat_example)
plot_regression_prediction(res$pred, res$response)
}

}
