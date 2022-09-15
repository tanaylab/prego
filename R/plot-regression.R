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
#' @description this would return a scatter plot of the response variable averages vs the regression model's prediction
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
    tibble(resp = response, pred = pred) %>%
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


#' Plot the comulative of the regression model's prediction stratified by the response variable
#'
#' @param pred the 'pred' field from the regression result
#' @param response the 'response' field from the regression result (the response variable). Should be binary (0/1).
#'
#' @examples
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], score_metric = "ks")
#' plot_regression_prediction_binary(res_binary$pred, res_binary$response)
#'
#' @export
plot_regression_prediction_binary <- function(pred, response) {
    if (!is_binary_response(response)) {
        cli_abort("Response is not binary")
    }

    # find the the threshold that maximizes the difference between the two cumulative distributions
    # in order to plot the KS D statistic
    pred_1 <- pred[response == 1]
    pred_0 <- pred[response == 0]
    cdf1 <- function(x) 1 - ecdf(pred_1)(x)
    cdf0 <- function(x) 1 - ecdf(pred_0)(x)
    min_max <- seq(min(c(pred_1, pred_0)), max(c(pred_1, pred_0)), length.out = length(pred))
    x0 <- min_max[which.max(abs(cdf0(min_max) - cdf1(min_max)))]
    y0 <- cdf0(x0)
    y1 <- cdf1(x0)
    ks <- suppressWarnings(ks.test(pred_1, pred_0))

    tibble(response = factor(response), pred = pred) %>%
        ggplot(aes(x = pred, y = 1 - ..y.., color = response)) +
        geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
            linetype = "dashed", color = "gray"
        ) +
        stat_ecdf() +
        scale_color_manual(name = "Response", values = c("1" = "red", "0" = "blue")) +
        theme_classic() +
        xlab("") +
        ylab("1 - ECDF") +
        labs(
            title = "Regression prediction",
            subtitle = as.expression(
                substitute(
                    italic("Kolmogorov-Smirnov D" ~ "=" ~ stat ~ "p-value" ~ "=" ~ pval),
                    list(
                        stat = paste0(round(ks$statistic, digits = 3), ","),
                        pval = round(ks$p.value, digits = 3)
                    )
                )
            )
        )
}

#' Plot the regression results
#'
#' @description
#'
#' Plot QC of the regression results
#'
#' @param reg output of \code{regress_pwm}
#' @param response the response variable
#' @param title a title for the plot (optional)
#'
#' @return a patchwork object
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_regression_qc(res)
#'
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' plot_regression_qc(res_binary)
#'
#' @export
plot_regression_qc <- function(reg, response = NULL, title = NULL) {
    if (is.null(response)) {
        if (!("response" %in% names(reg))) {
            cli_abort("{.field response} is missing from the regression result. Please provide one or set {.code include_response=TRUE} in {.code regress_pwm()}")
        }
        response <- reg$response
    }
    if (is_binary_response(response)) {
        p_pred <- plot_regression_prediction_binary(reg$pred, response)
    } else {
        p_pred <- plot_regression_prediction(reg$pred, response)
    }
    design <- "LS
               R#"
    p <- patchwork::wrap_plots(
        L = plot_pssm_logo(reg$pssm),
        S = plot_spat_model(reg$spat),
        R = p_pred,
        design = design
    )

    if (!is.null(title)) {
        p <- p + patchwork::plot_annotation(title = title)
    }
    return(p)
}
