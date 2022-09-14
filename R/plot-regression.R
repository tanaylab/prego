#' Plot LOGO of the pssm result from the regression
#'
#' @param pssm the 'pssm' field from the regression result
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_pssm_logo(res$pssm)
#'
#' @export
plot_pssm_logo <- function(pssm) {
    pfm <- pssm %>%
        as.data.frame() %>%
        tibble::column_to_rownames("pos") %>%
        as.matrix() %>%
        t()
    ggseqlogo::ggseqlogo(pfm) +
        ggtitle("Sequence model")
}

#' Plot spatial model of the regression result
#'
#' @param spat the 'spat' field from the regression result
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_spat_model(res$spat)
#'
#' @export
plot_spat_model <- function(spat) {
    spat %>%
        ggplot(aes(x = bin, y = spat_factor)) +
        geom_line() +
        geom_point() +
        theme_classic() +
        xlab("Position") +
        ylab("Spatial factor") +
        ggtitle("Spatial model")
}

#' Plot response variable averages vs the regression model's prediction
#'
#' @param pred the 'pred' field from the regression result
#' @param response the 'response' field from the regression result (the response variable)
#' @param point_size the size of the points in the plot (default: 0.5)
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_regression_prediction(res$pred, res$response)
#'
#' @export
plot_regression_prediction <- function(pred, response, point_size = 0.5) {
    if (is.matrix(response)) {
        response <- rowMeans(response)
    }
    r2 <- cor(response, pred)^2
    tibble::tibble(resp = response, pred = pred) %>%
        ggplot(aes(x = resp, y = pred)) +
        geom_point(size = point_size) +
        theme_classic() +
        xlab("Response") +
        ylab("Prediction") +
        labs(
            title = "Regression prediction",
            subtitle = as.expression(substitute(italic(r)^2 ~ "=" ~ r2, list(r2 = round(r2, digits = 3))))
        ) +
        theme(aspect.ratio = 1)
}

#' Plot the regression results
#'
#' @description
#'
#' Plot QC of the regression results
#'
#' @param reg output of \code{regress_pwm}
#'
#' @return a patchwork object
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_regression_qc(res)
#'
#' @export
plot_regression_qc <- function(reg) {
    design <- "LS
               R#"
    patchwork::wrap_plots(
        L = plot_pssm_logo(reg$pssm),
        S = plot_spat_model(reg$spat),
        R = plot_regression_prediction(reg$pred, reg$response),
        design = design
    )
}
