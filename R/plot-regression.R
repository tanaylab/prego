#' Plot spatial model of the regression result
#'
#' @param spat the 'spat' field from the regression result
#' @param title a title for the plot (optional)
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_spat_model(res$spat)
#' }
#'
#' @export
plot_spat_model <- function(spat, title = "Spatial model") {
    spat %>%
        ggplot(aes(x = bin, y = spat_factor)) +
        geom_line() +
        geom_point() +
        theme_classic() +
        xlab("Position") +
        ylab("Spatial factor") +
        ggtitle(title)
}

#' Plot response variable averages vs the regression model's prediction
#'
#' @description this would return a scatter plot of the response variable averages vs the regression model's prediction
#'
#' @param pred the 'pred' field from the regression result
#' @param response the 'response' field from the regression result (the response variable)
#' @param point_size the size of the points in the plot (default: 0.5)
#' @param alpha the transparency of the points in the plot (default: 1)
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_regression_prediction(res$pred, res$response)
#' }
#'
#' @export
plot_regression_prediction <- function(pred, response, point_size = 0.5, alpha = 1) {
    if (is.matrix(response)) {
        response <- rowMeans(response)
    }
    correlation <- cor(response, pred)
    r2 <- correlation^2

    tibble(resp = response, pred = pred) %>%
        ggplot(aes(x = resp, y = pred)) +
        geom_point(size = point_size, alpha = alpha) +
        theme_classic() +
        xlab("Response") +
        ylab("Prediction") +
        labs(
            title = "Regression prediction",
            subtitle = as.expression(substitute(italic(r)^2 ~ "=" ~ r2 ~ "," ~ italic(r) ~ "=" ~ correlation, list(r2 = round(r2, digits = 3), correlation = round(correlation, digits = 3))))
        ) +
        theme(aspect.ratio = 1)
}


#' Plot the cumulative of the regression model's prediction stratified by the response variable
#'
#' @param pred the 'pred' field from the regression result
#' @param response the 'response' field from the regression result (the response variable). Should be binary (0/1).
#'
#' @examples
#' \dontrun{
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], score_metric = "ks")
#' plot_regression_prediction_binary(res_binary$pred, res_binary$response)
#' }
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
    # Because we inverted the CDF, we need to flip the hypothesis test
    # (alternative = 'less' instead of 'greater')
    ks <- suppressWarnings(ks.test(pred_1, pred_0, alternative = "less"))

    tibble(response = factor(response), pred = pred) %>%
        ggplot(aes(x = pred, y = 1 - after_stat(y), color = response)) +
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
#' @param subtitle a subtitle for the plot (optional)
#' @param caption a caption for the plot (optional). When caption is NULL a default caption would be plotted.
#' @param point_size the size of the points in the scatter plot
#' @param alpha the transparency of the points in the scatter plot
#'
#'
#' @return a patchwork object
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_regression_qc(res)
#'
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], screen_db = TRUE)
#' plot_regression_qc(res_binary)
#' }
#'
#' @export
plot_regression_qc <- function(reg,
                               response = NULL,
                               title = glue("Motif regression results (consensus: {reg$consensus})"),
                               subtitle = NULL,
                               caption = NULL,
                               point_size = 0.5,
                               alpha = 0.5) {
    if (is.null(response)) {
        if (!("response" %in% names(reg))) {
            cli_abort("{.field response} is missing from the regression result. Please provide one or set {.code include_response=TRUE} in {.code regress_pwm()}")
        }
        response <- reg$response
    }
    if (is.null(caption)) {
        caption <- glue("# of 1: {sum(response == 1)}, \\
                        # of 0: {sum(response == 0)}, \\
                        seed: {reg$seed_motif}")
    }
    if (is_binary_response(response)) {
        p_pred <- plot_regression_prediction_binary(reg$pred, response)
    } else {
        p_pred <- plot_regression_prediction(reg$pred, response, point_size = point_size, alpha = alpha)
    }

    if (is.null(reg$db_match) || is.null(reg$db_match_cor)) {
        match_df <- pssm_match(reg$pssm, all_motif_datasets())
        best_motif <- match_df[1, ]
        reg$db_match <- best_motif$motif
        reg$db_match_cor <- best_motif$cor
    }

    m_subtitle <- glue("cor: {round(reg$db_match_cor, digits = 3)}")

    if (is_binary_response(response)) {
        if (!is.null(reg$db_match_ks)) {
            m_subtitle <- glue("cor: {round(reg$db_match_cor, digits = 3)}, KS D: {round(reg$db_match_ks$statistic, digits = 3)}")
        }
    } else {
        if (!is.null(reg$db_match_r2)) {
            m_subtitle <- glue("cor: {round(reg$db_match_cor, digits = 3)}, R2: {round(reg$db_match_r2, digits = 3)}")
        }
    }

    p_match <- plot_pssm_logo_dataset(
        reg$db_match,
        all_motif_datasets(),
        title = glue("Most similar to:\n{reg$db_match}"),
        subtitle = m_subtitle
    )

    if (!is.null(reg$db_motif)) {
        p_db_logo <- plot_pssm_logo_dataset(
            reg$db_motif,
            all_motif_datasets(),
            title = glue("Best motif in database:\n{reg$db_motif}"),
            subtitle = glue("score: {round(reg$db_motif_score, digits = 3)}")
        )

        design <- "LS
                   MR
                   DN"
        p <- patchwork::wrap_plots(
            L = plot_pssm_logo(reg$pssm),
            S = plot_spat_model(reg$spat),
            R = p_pred,
            M = p_match,
            D = p_db_logo,
            design = design
        )
    } else {
        design <- "LS
               MR"
        p <- patchwork::wrap_plots(
            L = plot_pssm_logo(reg$pssm),
            S = plot_spat_model(reg$spat),
            R = p_pred,
            M = p_match,
            design = design
        )
    }


    if (!is.null(title) || !is.null(subtitle) || !is.null(caption)) {
        p <- p + patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption)
    }
    return(p)
}

#' Plot the regression results for multiple motifs
#'
#' @description plot the regression results when \code{motif_num} > 1
#'
#' @examples
#' \dontrun{
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 3], motif_num = 3)
#' plot_regression_qc_multi(res_binary)
#' }
#'
#' @inheritParams plot_regression_qc
#' @export
plot_regression_qc_multi <- function(reg, title = glue("Motif regression results (consensus: {reg$consensus})"),
                                     subtitle = NULL,
                                     caption = NULL,
                                     point_size = 0.01,
                                     alpha = 0.5,
                                     response = NULL) {
    if (!("models" %in% names(reg)) || !("multi_stats" %in% names(reg))) {
        cli_abort("The regression result does not contain multiple motifs")
    }

    if (is.null(response)) {
        if (!is.null(reg$response)) {
            response <- reg$response
        } else if (!is.null(reg$models[[1]]$response)) {
            response <- reg$models[[1]]$response
        } else {
            cli_abort("{.field response} is missing from the regression result. Please provide one or set {.code include_response=TRUE} in {.code regress_pwm()}")
        }
    }


    if (is.null(caption)) {
        caption <- glue("# of 1: {sum(response == 1)}, \\
                        # of 0: {sum(response == 0)}, \\
                        seed: {reg$seed_motif}")
    }

    models <- reg$models

    spatial_p <- purrr::imap(models, ~ plot_spat_model(.x$spat))
    motifs_p <- purrr::imap(models, ~ plot_pssm_logo(.x$pssm, title = paste0("Motif #", .y), subtitle = glue("score = {round(reg$multi_stats$score[.y], digits=3)}, combined score = {round(reg$multi_stats$comb_score[.y], digits=3)}")))
    motifs_match_p <- purrr::imap(models, ~ {
        if (is.null(.x$db_match_pssm)) {
            return(NULL)
        }
        if (!is.null(.x$db_match_ks)) {
            plot_pssm_logo(.x$db_match_pssm, title = .x$db_match, subtitle = glue("KS D: {round(.x$db_match_ks$statistic, digits = 3)}"))
        } else {
            plot_pssm_logo(.x$db_match_pssm, title = .x$db_match, subtitle = glue("R^2: {round(.x$db_match_r2, digits = 3)}"))
        }
    })

    prediction_p <- purrr::imap(models, ~ {
        if (is_binary_response(response)) {
            plot_regression_prediction_binary(.x$pred, response)
        } else {
            plot_regression_prediction(.x$pred, response, point_size = point_size, alpha = alpha)
        }
    })


    p_scores <- reg$multi_stats %>%
        dplyr::rename(`Score` = score, `Combined score` = comb_score) %>%
        tidyr::pivot_longer(cols = c("Score", "Combined score")) %>%
        mutate(label = round(value, digits = 3)) %>%
        ggplot(aes(x = factor(model), y = value, label = label)) +
        geom_col() +
        theme_classic() +
        geom_text(vjust = -0.5) +
        scale_y_continuous(expand = expansion(mult = 0.1)) +
        xlab("Model") +
        ylab("Score") +
        facet_wrap(. ~ name, scales = "free_y", ncol = 2)

    p <- patchwork::wrap_plots(
        A = patchwork::wrap_plots(spatial_p, ncol = 1),
        B = patchwork::wrap_plots(motifs_p, ncol = 1),
        C = patchwork::wrap_plots(motifs_match_p, ncol = 1),
        D = patchwork::wrap_plots(prediction_p, ncol = 1),
        E = p_scores,
        design = "ABCD
                  EEEE",
        heights = c(0.8, 0.2)
    )

    if (!is.null(title) || !is.null(subtitle) || !is.null(caption)) {
        p <- p + patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption)
    }

    return(p)
}
