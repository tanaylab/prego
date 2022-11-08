#' Run PWM regression on a sample of the data
#'
#' @description The optimization would be performed with a sampled dataset of size \code{sample_frac}, or explicit sampled indices \code{sample_idxs}. Note that \code{multi_kmers} is TRUE by default.
#'
#' @param sample_frac fraction of the dataset to sample. When \code{response} is categorical (0 and 1), the sampling would be stratified by the category, i.e. \code{sample_frac} can be a vector of length 2 with the fraction of 0 and 1 responses to sample respectively.
#' If NULL - the default would be 0.1 for continuous variables, and for binary variables - the number of 0 responses would be equal to \code{sample_ratio} times the number of 1 responses.
#' @param sample_idxs indices of the sequences to use. If NULL, the indices would be sampled using \code{sample_frac}.
#' @param sample_ratio ratio between the '1' category and the '0' category in the sampled dataset. Relevant only when \code{sample_frac} is NULL.
#'
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm.sample(cluster_sequences_example, cluster_mat_example[, 1], final_metric = "ks", screen_db = TRUE)
#' res$pssm
#' res$spat
#' head(res$pred)
#'
#' plot_regression_qc(res)
#' }
#'
#' @inheritParams regress_pwm
#' @inheritDotParams regress_pwm
#' @inheritDotParams screen_kmers
#' @inherit regress_pwm return
#' @export
regress_pwm.sample <- function(sequences,
                               response,
                               bidirect = TRUE,
                               spat_min = 0,
                               spat_max = NULL,
                               include_response = TRUE,
                               motif_num = 1,
                               multi_kmers = TRUE,
                               sample_frac = NULL,
                               sample_idxs = NULL,
                               sample_ratio = 1,
                               parallel = getOption("prego.parallel", TRUE),
                               match_with_db = FALSE,
                               screen_db = FALSE,
                               motif_dataset = all_motif_datasets(),
                               seed = 60427,
                               final_metric = NULL,
                               unif_prior = 0.05,
                               ...) {
    set.seed(seed)
    if (is.null(nrow(response))) {
        response <- matrix(response, ncol = 1)
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    cli_alert_info("Performing sampled optimization")
    if (is.null(sample_idxs)) {
        sample_idxs <- sample_response(response, sample_frac, sample_ratio, seed)
    }

    sequences_s <- sequences[sample_idxs]
    response_s <- response[sample_idxs, , drop = FALSE]

    res <- regress_pwm(
        sequences = sequences_s,
        response = response_s,
        bidirect = bidirect,
        unif_prior = unif_prior,
        spat_min = spat_min,
        spat_max = spat_max,
        motif_num = motif_num,
        multi_kmers = multi_kmers,
        include_response = FALSE,
        verbose = FALSE,
        match_with_db = FALSE,
        parallel = parallel,
        final_metric = final_metric,
        seed = seed,
        ...
    )

    res$sample_idxs <- sample_idxs

    # fill predictions for all the sequences
    res$pred <- compute_pwm(sequences, res$pssm, res$spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect)

    if (motif_num > 1 && "models" %in% names(res)) {
        res$models <- purrr::map(res$models, ~ {
            .x$pred <- compute_pwm(sequences, .x$pssm, .x$spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect)
            .x
        })
    }

    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])], alternative = "less"))
    }

    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset, parallel = parallel)
    }

    if (screen_db) {
        res <- add_regression_db_screen(res, response, sequences, motif_dataset, final_metric, prior = unif_prior, bidirect = bidirect, parallel = parallel)
    }

    cli_alert_success("Finished running regression. Consensus: {.val {res$consensus}}")

    if (is_binary_response(response)) {
        cli_alert_success("KS test D: {.val {round(res$ks$statistic, digits=4)}}, p-value: {.val {res$ks$p.value}}")
    } else {
        cli_alert_success("R^2: {.val {round(res$r2, digits=4)}}")
    }

    return(res)
}

sample_response <- function(response, sample_frac = NULL, sample_ratio = 1, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (is.null(sample_frac)) {
        if (is_binary_response(response)) {
            sample_frac <- c(pmin(1, sample_ratio * sum(response[, 1] == 1) / sum(response[, 1] == 0)), 1)
        } else {
            cli_alert_info("Using {.code sample_frac = 0.1}")
            sample_frac <- 0.1
        }
    }
    cli_alert_info("Sampling {.val {round(sample_frac, digits = 2)}} of the dataset")
    categorical <- ncol(response) == 1 && all(response[, 1] == 0 | response[, 1] == 1)
    if (categorical) {
        cli_alert_info("Stratified sampling")
        if (length(sample_frac) == 1) {
            sample_frac <- c(sample_frac, sample_frac)
        }
        samp_idx_0 <- sample(which(response[, 1] == 0), size = round(sample_frac[1] * sum(response[, 1] == 0)))
        samp_idx_1 <- sample(which(response[, 1] == 1), size = round(sample_frac[2] * sum(response[, 1] == 1)))
        sample_idxs <- c(samp_idx_0, samp_idx_1)
    } else {
        sample_idxs <- sample(1:nrow(response), size = floor(sample_frac * nrow(response)))
    }

    if (is_binary_response(response)) {
        cli_alert_info("sampled {.val {sum(response[sample_idxs, 1] == 0)}} zeros and {.val {sum(response[sample_idxs, 1] == 1)}} ones")
    }
    return(sample_idxs)
}
