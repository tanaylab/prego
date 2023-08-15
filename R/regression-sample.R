#' Run PWM regression on a sample of the data
#'
#' @description The optimization would be performed with a sampled dataset of size \code{sample_frac}, or explicit sampled indices \code{sample_idxs}.
#'
#' @param sample_frac fraction of the dataset to sample. When \code{response} is categorical (0 and 1), the sampling would be stratified by the category, i.e. \code{sample_frac} can be a vector of length 2 with the fraction of 0 and 1 responses to sample respectively.
#' If NULL - the default would be 0.1 for continuous variables, and for binary variables - the number of 0 responses would be equal to \code{sample_ratio} times the number of 1 responses.
#' @param sample_idxs indices of the sequences to use. If NULL, the indices would be sampled using \code{sample_frac}.
#' @param sample_ratio ratio between the '1' category and the '0' category in the sampled dataset. Relevant only when \code{sample_frac} is NULL.
#'
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm.sample(
#'     cluster_sequences_example,
#'     cluster_mat_example[, 1],
#'     final_metric = "ks",
#'     screen_db = TRUE
#' )
#'
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
                               spat_bin_size = 40,
                               spat_num_bins = 7,
                               bidirect = TRUE,
                               include_response = TRUE,
                               motif_num = 1,
                               multi_kmers = TRUE,
                               sample_frac = NULL,
                               sample_idxs = NULL,
                               sample_ratio = 1,
                               parallel = getOption("prego.parallel", TRUE),
                               match_with_db = TRUE,
                               screen_db = FALSE,
                               motif_dataset = all_motif_datasets(),
                               seed = 60427,
                               final_metric = NULL,
                               unif_prior = 0.05,
                               alternative = "two.sided",
                               energy_func = NULL,
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
        spat_bin_size = spat_bin_size,
        spat_num_bins = spat_num_bins,
        motif_num = motif_num,
        multi_kmers = multi_kmers,
        include_response = FALSE,
        verbose = FALSE,
        match_with_db = FALSE,
        screen_db = FALSE,
        parallel = parallel,
        final_metric = final_metric,
        seed = seed,
        alternative = alternative,
        sample_for_kmers = FALSE,
        energy_func = energy_func,
        ...
    )

    res$sample_idxs <- sample_idxs

    spat <- calc_spat_min_max(spat_bin_size, spat_num_bins, nchar(sequences_s[1]))

    # fill predictions for all the sequences
    res$pred <- compute_pwm(sequences, res$pssm, res$spat, spat_min = spat$spat_min, spat_max = spat$spat_max, bidirect = bidirect)

    if (motif_num > 1 && "models" %in% names(res)) {
        res$models <- purrr::map(res$models, ~ {
            .x$pred <- compute_pwm(sequences, .x$pssm, .x$spat, spat_min = spat$spat_min, spat_max = spat$spat_max, bidirect = bidirect)
            .x
        })
    }

    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])], alternative = alternative))
    }

    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset, parallel = parallel, alternative = alternative)
    }

    if (screen_db) {
        res <- add_regression_db_screen(res, response, sequences, motif_dataset, final_metric, prior = unif_prior, bidirect = bidirect, parallel = parallel, alternative = alternative)
    }

    cli_alert_success("Finished running regression. Consensus: {.val {res$consensus}}")

    if (is_binary_response(response)) {
        cli_alert_success("KS test D: {.val {round(res$ks$statistic, digits=4)}}, p-value: {.val {res$ks$p.value}}")
    } else {
        cli_alert_success("R^2: {.val {round(res$r2, digits=4)}}")
    }

    res <- add_predict_function(res, spat, bidirect, energy_func)

    res$spat_min <- spat$spat_min
    res$spat_max <- spat$spat_max
    res$spat_bin_size <- spat_bin_size
    res$bidirect <- bidirect
    res$seq_length <- nchar(sequences[1])

    return(res)
}
