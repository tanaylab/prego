#' Cross-validate a two-phase regression model
#'
#' @description Perform cross-validation on a two-phase regression model. You can either provide explicit folds, or use the \code{nfolds} argument to set the number of folds. If the response is binary (0 and 1) or a \code{categories} vector is given, the folds would be stratified by the response/categories.
#'
#' @param nfolds number of folds for cross-validation. Can be NULL if \code{folds} are provided.
#' @param metric metric to use for cross-validation. One of 'ks' or 'r2'. If NULL - 'ks' would be set for binary response and 'r2' for continuous response.
#' @param folds vector of fold numbers for each sequence (optional)
#' @param categories vector of categories for each sequence (optional)
#' @param parallel whether to run the cross-validation in parallel.
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{cv_models: }{a list of models, one for each fold.}
#' \item{cv_pred: }{a vector of predictions for each sequence.}
#' \item{score: }{score of the model on the cross-validated predictions.}
#' \item{cv_scores: }{a vector of scores for each fold.}
#' \item{folds: }{a vector with the fold number for each sequence.}
#' }
#'
#' @inheritParams regress_pwm
#'
#' @export
regress_pwm_two_phase.cv <- function(sequences,
                                     response,
                                     nfolds = NULL,
                                     metric = NULL,
                                     folds = NULL,
                                     categories = NULL,
                                     seed = 60427,
                                     parallel = getOption("prego.parallel"),
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

    if ("first_phase_idxs" %in% names(list(...))) {
        cli_abort("The {.field first_phase_idxs} argument is not supported in {.fun regress_pwm_two_phase.cv}")
    }

    cv_res <- plyr::llply(unique(folds), function(f) {
        cli_h1("Cross-validation fold {.val {f}}")
        train_idxs <- which(folds != f)
        test_idxs <- which(folds == f)
        res <- regress_pwm_two_phase(sequences[train_idxs],
            response[train_idxs, , drop = FALSE],
            seed = seed,
            ...
        )
        res$test_pred <- compute_pwm(sequences[test_idxs], res$pssm, res$spat)
        res$test_score <- score_model(metric, response[test_idxs, , drop = FALSE], res$test_pred)
        return(res)
    }, .parallel = parallel)

    names(cv_res) <- paste0("fold", unique(folds))
    cv_pred <- purrr::map_dfr(unique(folds), ~
        tibble(idx = which(folds == .x), pred = cv_res[[paste0("fold", .x)]]$test_pred)) %>%
        arrange(idx) %>%
        pull(pred)
    score <- score_model(metric, response, cv_pred)

    cli_alert_success("Cross-validation score: {.val {score}}")

    res <- list(
        cv_models = cv_res,
        cv_pred = cv_pred,
        score = score,
        cv_scores = sapply(cv_res, function(x) x$test_score),
        folds = folds
    )

    return(res)
}

score_model <- function(metric, response, pred) {
    if (metric == "ks") {
        score <- suppressWarnings(ks.test(pred[response == 0], pred[response == 1])$statistic)
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
