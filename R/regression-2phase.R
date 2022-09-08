#' Run a 2-phase PWM regression
#'
#' @description The first phase of the optimization would be performed with a sampled dataset of size \code{two_phase_sample_frac} and then the optimization would be performed on the full dataset while initializing the motif from the sampled dataset. You can also give explicit indices of the sequences to use in the first phase using \code{first_phase_idxs}.
#'
#' @param two_phase_sample_frac fraction of the dataset to sample for the first phase of the optimization (default: 0.1). When \code{response} is categorical (0 and 1), the sampling would be stratified by the category,
#' i.e. \code{two_phase_sample_frac} can be a vector of length 2 with the fraction of 0 and 1 responses to sample
#' respectively.
#' @param first_phase_idxs indices of the sequences to use in the first phase of the optimization. If NULL, the indices would be sampled using \code{two_phase_sample_frac}.
#'
#' @examples
#' res <- regress_pwm_two_phase(sequences_example, response_mat_example, two_phase_sample_frac = 0.1)
#' res$pssm
#' res$spat
#' head(res$pred)
#'
#' @inheritParams regress_pwm
#' @export
regress_pwm_two_phase <- function(sequences,
                                  response,
                                  motif = NULL,
                                  motif_length = 15,
                                  score_metric = "r2",
                                  bidirect = TRUE,
                                  spat_min = 0,
                                  spat_max = NULL,
                                  spat_bin = 50,
                                  improve_epsilon = 0.0001,
                                  min_nuc_prob = 0.001,
                                  unif_prior = 0.05,
                                  is_train = NULL,
                                  include_response = TRUE,
                                  seed = 60427,
                                  verbose = FALSE,
                                  kmer_length = 8,
                                  two_phase_sample_frac = 0.1,
                                  first_phase_idxs = NULL,
                                  ...) {
    if (is.null(nrow(response))) {
        response <- matrix(response, ncol = 1)
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    cli_alert_info("Performing two phase optimization")
    if (is.null(first_phase_idxs)) {
        cli_alert_info("Sampling {.val {two_phase_sample_frac}} of the dataset for the first phase")
        categorical <- ncol(response) == 1 && all(response[, 1] == 0 | response[, 1] == 1)
        if (categorical) {
            cli_alert_info("Stratified sampling")
            if (length(two_phase_sample_frac) == 1) {
                two_phase_sample_frac <- c(two_phase_sample_frac, two_phase_sample_frac)
            }
            samp_idx_0 <- sample(which(response[, 1] == 0), size = round(two_phase_sample_frac[1] * sum(response[, 1] == 0)))
            samp_idx_1 <- sample(which(response[, 1] == 1), size = round(two_phase_sample_frac[2] * sum(response[, 1] == 1)))
            first_phase_idxs <- c(samp_idx_0, samp_idx_1)
        } else {
            first_phase_idxs <- sample(1:length(sequences), size = floor(two_phase_sample_frac * length(sequences)))
        }
    }

    sequences_s <- sequences[first_phase_idxs]
    response_s <- response[first_phase_idxs, , drop = FALSE]
    cli_alert_info("sampled {.val {sum(response_s[, 1] == 0)}} 0s and {.val {sum(response_s[, 1] == 1)}} 1s")

    cli_alert_info("Running regression on the sampled dataset")
    res_s <- regress_pwm(
        sequences = sequences_s,
        response = response_s,
        motif = NULL,
        motif_length = motif_length,
        score_metric = score_metric,
        bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_bin = spat_bin,
        improve_epsilon = improve_epsilon,
        min_nuc_prob = min_nuc_prob,
        unif_prior = unif_prior,
        is_train = is_train,
        include_response = include_response,
        seed = seed,
        verbose = verbose,
        kmer_length = kmer_length,
        motif_num = 1
    )

    cli_alert_info("Running regression on the full dataset")
    res <- regress_pwm(
        sequences = sequences,
        response = response,
        motif = res_s$pssm,
        motif_length = motif_length,
        score_metric = score_metric,
        bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_bin = spat_bin,
        improve_epsilon = improve_epsilon,
        min_nuc_prob = min_nuc_prob,
        unif_prior = unif_prior,
        is_train = is_train,
        include_response = include_response,
        seed = seed,
        verbose = verbose,
        kmer_length = kmer_length,
        motif_num = 1
    )

    return(res)
}
