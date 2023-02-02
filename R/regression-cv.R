#' Cross-validate a PWM regression model
#'
#' @description Perform cross-validation on a PWM regression model. You can either provide explicit folds, or use the \code{nfolds} argument to set the number of folds. If the response is binary (0 and 1) or a \code{categories} vector is given, the folds would be stratified by the response/categories.
#'
#' @param nfolds number of folds for cross-validation. Can be NULL if \code{folds} are provided.
#' @param metric metric to use for cross-validation. One of 'ks' or 'r2'. If NULL - 'ks' would be set for binary response and 'r2' for continuous response.
#' @param folds vector of fold numbers for each sequence (optional)
#' @param categories vector of categories for each sequence (optional)
#' @param use_sample whether to use sampled optimization or not.
#' @param parallel whether to run the cross-validation in parallel.
#' @param fold_parallel whether to run the optimization in each fold in parallel. It is recommended to set this to FALSE if \code{parallel} is TRUE.
#' @param add_full_model whether to add the full model (without cross-validation) to the results.
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{cv_models: }{a list of models, one for each fold.}
#' \item{cv_pred: }{a vector of predictions for each sequence.}
#' \item{score: }{score of the model on the cross-validated predictions.}
#' \item{cv_scores: }{a vector of scores for each fold.}
#' \item{folds: }{a vector with the fold number for each sequence.}
#' \item{full_model: }{The full model (without cross-validation), if \code{add_full_model} is TRUE.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm.cv(
#'     cluster_sequences_example, cluster_mat_example[, 1],
#'     nfolds = 5, use_sample = TRUE, sample_frac = c(0.1, 1)
#' )
#' res$score
#' res$cv_scores
#'
#' plot(
#'     res$cv_pred,
#'     res$full_model$pred,
#'     xlab = "CV predictions", ylab = "Full model predictions", cex = 0.1
#' )
#' plot_regression_prediction_binary(res$cv_pred, cluster_mat_example[, 1])
#' plot_regression_prediction_binary(res$full_model$pred, cluster_mat_example[, 1])
#'
#' # without sampling
#' res <- regress_pwm.cv(
#'     cluster_sequences_example, cluster_mat_example[, 1],
#'     nfolds = 5, use_sample = FALSE
#' )
#' res$score
#' res$cv_scores
#' plot(res$cv_pred,
#'     res$full_model$pred,
#'     xlab = "CV predictions", ylab = "Full model predictions", cex = 0.1
#' )
#' }
#'
#' @inheritParams regress_pwm
#' @inheritDotParams regress_pwm
#' @inheritDotParams regress_pwm.sample
#'
#' @export
regress_pwm.cv <- function(sequences,
                           response,
                           nfolds = NULL,
                           metric = NULL,
                           folds = NULL,
                           categories = NULL,
                           use_sample = FALSE,
                           seed = 60427,
                           parallel = getOption("prego.parallel", FALSE),
                           fold_parallel = !parallel && getOption("prego.parallel", FALSE),
                           add_full_model = TRUE,
                           alternative = "less",
                           ...) {
    set.seed(seed)
    if (is.null(folds)) {
        folds <- get_cv_folds(response, nfolds, categories)
    } else {
        if (length(folds) != length(sequences)) {
            cli_abort("The number of folds and the number of sequences do not match")
        }
    }

    if (!is.matrix(response)) {
        response <- matrix(response, ncol = 1)
    }

    if (is.null(metric)) {
        if (is_binary_response(response)) {
            cli_alert_info("Response is binary: setting metric to {.val ks}")
            metric <- "ks"
        } else {
            cli_alert_info("Response is continuous: setting metric to {.val r2}")
            metric <- "r2"
        }
    }

    if (use_sample) {
        cli_alert_info("Using sampled optimization")
        func <- regress_pwm.sample
        if ("sample_idxs" %in% names(list(...))) {
            cli_abort("The {.field sample_idxs} argument is not supported in {.fun regress_pwm.cv}")
        }
    } else {
        func <- regress_pwm
    }

    regression_func <- purrr::partial(func, alternative = alternative, parallel = fold_parallel)

    cv_res <- plyr::llply(unique(folds), function(f) {
        cli_h1("Cross-validation fold {.val {f}}")
        train_idxs <- which(folds != f)
        test_idxs <- which(folds == f)
        res <- regression_func(sequences[train_idxs],
            response[train_idxs, , drop = FALSE],
            seed = seed,
            ...
        )

        res$test_pred <- res$predict(sequences[test_idxs])

        res$test_score <- score_model(metric, response[test_idxs, , drop = FALSE], res$test_pred, alternative = alternative)
        return(res)
    }, .parallel = parallel)

    names(cv_res) <- paste0("fold", unique(folds))
    cv_pred <- purrr::map_dfr(unique(folds), ~
        tibble(idx = which(folds == .x), pred = cv_res[[paste0("fold", .x)]]$test_pred)) %>%
        arrange(idx) %>%
        pull(pred)
    score <- score_model(metric, response, cv_pred, alternative = alternative)

    cli_alert_success("Cross-validation score: {.val {score}}")

    res <- list(
        cv_models = cv_res,
        cv_pred = cv_pred,
        score = score,
        cv_scores = sapply(cv_res, function(x) x$test_score),
        folds = folds
    )

    if (add_full_model) {
        res$full_model <- regression_func(sequences,
            response,
            seed = seed,            
            ...
        )
    }

    return(res)
}

score_model <- function(metric, response, pred, alternative = "two.sided") {
    if (metric == "ks") {
        score <- suppressWarnings(ks.test(pred[response == 1], pred[response == 0], alternative = alternative)$statistic)
    } else if (metric == "r2") {
        score <- cor(response, pred)^2
    } else {
        cli_abort("Unknown metric {.val {metric}}")
    }
    return(score)
}


get_cv_folds <- function(response, nfolds, categories = NULL) {
    if (is.matrix(response)) {
        response <- response[, 1]
    }

    if (is.null(nfolds)) {
        cli_abort("Either {.field nfolds} or {.field folds} must be provided")
    }

    if (is_binary_response(response) || !is.null(categories)) {
        if (is.null(categories)) {
            if (is.matrix(response)) {
                categories <- response[, 1]
            } else {
                categories <- response
            }
        }
        cli_alert_info("Stratified sampling")
        folds <- tibble::tibble(idx = seq_along(response), cat = categories) %>%
            slice(sample(idx)) %>%
            group_by(cat) %>%
            mutate(fold = dplyr::ntile(idx, n = nfolds)) %>%
            ungroup() %>%
            arrange(idx) %>%
            pull(fold)
    } else {
        folds <- dplyr::ntile(seq_along(response), nfolds)
    }

    return(folds)
}
