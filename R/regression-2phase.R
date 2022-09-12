#' Run a 2-phase PWM regression
#'
#' @description The first phase of the optimization would be performed with a sampled dataset of size \code{two_phase_sample_frac} where different candidates of kmers would be regressed in order to find the best seed. Thenm the optimization would be performed on the full dataset while initializing the motif from the sampled dataset. You can also give explicit indices of the sequences to use in the first phase using \code{first_phase_idxs}.
#'
#' @param two_phase_sample_frac fraction of the dataset to sample for the first phase of the optimization (default: 0.1). When \code{response} is categorical (0 and 1), the sampling would be stratified by the category,
#' i.e. \code{two_phase_sample_frac} can be a vector of length 2 with the fraction of 0 and 1 responses to sample
#' respectively.
#' @param first_phase_idxs indices of the sequences to use in the first phase of the optimization. If NULL, the indices would be sampled using \code{two_phase_sample_frac}.
#' @param first_phase_metric metric to use in order to choose the best motif in the first phase of the optimization. One of 'ks' or 'r2'. Note that unlike \code{score_metric} which is used in the regression itself, this metric is used only for choosing the best motif in the first phase of the optimization out of all the runs on the sampled dataset.
#' @param kmer_length a vector of kmer lengths to screen in order to find the best seed motif.
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
                                  kmer_length = 5:8,
                                  min_gap = 0,
                                  max_gap = 1,
                                  min_kmer_cor = 0.1,
                                  two_phase_sample_frac = 0.1,
                                  first_phase_idxs = NULL,
                                  first_phase_metric = "r2",
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

    cli_h1("First phase")
    cli_h2("Generate candidate kmers")
    cand_kmers <- get_cand_kmers(sequences_s, response_s, kmer_length, min_gap, max_gap, min_kmer_cor, ...)

    cli_h2("Regress each candidate kmer on sampled data")
    cli_alert_info("Running {.val {length(cand_kmers)}} candidate kmers")
    res_s_list <- lapply(cand_kmers, function(motif) {
        regress_pwm(sequences_s,
            response_s,
            motif = motif,
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
            verbose = verbose
        )
    })
    scores <- sapply(res_s_list, function(x) x[[first_phase_metric]]$statistic)
    res_s <- res_s_list[[which.max(scores)]]

    cli_alert_info("Best motif in the first phase: {.val {res_s$seed_motif}}")

    cli_h1("Running regression on the full dataset")
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

    res$seed_motif <- res_s$seed_motif
    res$kmers <- res_s$kmers

    return(res)
}

get_cand_kmers <- function(sequences, response, kmer_length, min_gap, max_gap, min_kmer_cor, ...) {
    params <- expand.grid(kmer_length, min_gap:max_gap)
    colnames(params) <- c("len", "gap")
    all_kmers <- purrr::map_dfr(1:nrow(params), function(i) {
        screen_kmers(sequences, response, kmer_length = params$len[i], min_gap = 0, max_gap = params$gap[i], ...) %>%
            mutate(len = params$len[i], gap = params$gap[i])
    })

    best_motif <- all_kmers$kmer[which.max(abs(all_kmers$max_r2))] # return at least one motif

    cands <- all_kmers %>%
        filter(sqrt(max_r2) > min_kmer_cor) %>%
        group_by(len, gap) %>%
        slice_max(n = 2, order_by = abs(max_r2)) %>%
        pull(kmer)

    cands <- unique(c(best_motif, cands))

    return(cands)
}
