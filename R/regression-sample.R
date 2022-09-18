#' Run PWM regression on a sample of the data
#'
#' @description The optimization would be performed with a sampled dataset of size \code{sample_frac} (or explicit sampled indices \code{sample_idxs}) where different candidates of kmers would be regressed in order to find the best seed according to \code{final_metric}.
#'
#' @param sample_frac fraction of the dataset to sample. When \code{response} is categorical (0 and 1), the sampling would be stratified by the category, i.e. \code{sample_frac} can be a vector of length 2 with the fraction of 0 and 1 responses to sample respectively.
#' If NULL - the default would be 0.1 for continuous variables, and for binary variables - the number of 0 responses would be equal to \code{sample_ratio} times the number of 1 responses.
#' @param sample_idxs indices of the sequences to use. If NULL, the indices would be sampled using \code{sample_frac}.
#' @param sample_ratio ratio between the '1' category and the '0' category in the sampled dataset. Relevant only when \code{sample_frac} is NULL.
#' @param final_metric metric to use in order to choose the best motif. One of 'ks' or 'r2'. Note that unlike \code{score_metric} which is used in the regression itself, this metric is used only for choosing the best motif out of all the runs on the sampled dataset.
#' @param kmer_length a vector of kmer lengths to screen in order to find the best seed motif.
#' @param max_cands maximum number of kmer candidates to try.
#' @param verbose verbosity of the optimization.
#' @param parallel whether to run optimization in parallel. use \code{set_parallel}
#' to set the number of cores to use.
#'
#'
#' @examples
#' res <- regress_pwm.sample(cluster_sequences_example, cluster_mat_example[, 1], final_metric = "ks")
#' res$pssm
#' res$spat
#' head(res$pred)
#'
#' plot_regression_qc(res)
#'
#' @inheritParams regress_pwm
#' @inheritDotParams screen_kmers
#' @inherit regress_pwm return
#' @export
regress_pwm.sample <- function(sequences,
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
                               kmer_length = 6:8,
                               max_cands = 10,
                               min_gap = 0,
                               max_gap = 1,
                               min_kmer_cor = 0.1,
                               consensus_single_thresh = 0.6,
                               consensus_double_thresh = 0.85,
                               sample_frac = NULL,
                               sample_idxs = NULL,
                               sample_ratio = 1,
                               final_metric = "r2",
                               parallel = getOption("prego.parallel", FALSE),
                               match_with_db = FALSE,
                               motif_dataset = all_motif_datasets(),
                               ...) {
    set.seed(seed)
    if (is.null(nrow(response))) {
        response <- matrix(response, ncol = 1)
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (is.null(sample_frac)) {
        if (is_binary_response(response)) {
            sample_frac <- c(pmin(1, sample_ratio * sum(response[, 1] == 1) / sum(response[, 1] == 0)), 1)
        } else {
            cli_alert_info("Using {.code sample_frac = 0.1}")
            sample_frac <- 0.1
        }
    }

    cli_alert_info("Performing sampled optimization")

    if (is.null(sample_idxs)) {
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
            sample_idxs <- sample(seq_along(sequences), size = floor(sample_frac * length(sequences)))
        }
    }

    sequences_s <- sequences[sample_idxs]
    response_s <- response[sample_idxs, , drop = FALSE]
    if (is_binary_response(response_s)) {
        cli_alert_info("sampled {.val {sum(response_s[, 1] == 0)}} zeros and {.val {sum(response_s[, 1] == 1)}} ones")
    }

    cli_h3("Generate candidate kmers")
    cand_kmers <- get_cand_kmers(sequences_s, response_s, kmer_length, min_gap, max_gap, min_kmer_cor, verbose, parallel, max_cands = max_cands, ...)

    cli_h3("Regress each candidate kmer on sampled data")
    cli_alert_info("Running regression on {.val {length(cand_kmers)}} candidate kmers")
    cli_ul(c(
        "Bidirectional: {.val {bidirect}}",
        "Spat min: {.val {spat_min}}",
        "Spat max: {.val {spat_max}}",
        "Spat bin: {.val {spat_bin}}",
        "Improve epsilon: {.val {improve_epsilon}}",
        "Min nuc prob: {.val {min_nuc_prob}}",
        "Uniform prior: {.val {unif_prior}}",
        "Score metric: {.val {score_metric}}",
        "Seed: {.val {seed}}"
    ))
    res_s_list <- plyr::llply(cli_progress_along(cand_kmers), function(i) {
        motif <- cand_kmers[i]
        cli_alert("regressing with seed: {.val {motif}}")
        r <- regress_pwm(sequences_s,
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
            include_response = FALSE,
            seed = seed,
            verbose = FALSE,
            consensus_single_thresh = consensus_single_thresh,
            consensus_double_thresh = consensus_double_thresh,
            match_with_db = FALSE
        ) %>%
            suppressMessages()
        if (final_metric == "ks") {
            if (!is_binary_response(response_s)) {
                cli_abort("Cannot use {.field final_metric} {.val ks} when {.field response} is not binary")
            }
            r$score <- r[[final_metric]]$statistic
        } else if (final_metric == "r2") {
            r$score <- r[[final_metric]]
        } else {
            cli_abort("Unknown {.field final_metric} (can be 'ks' or 'r2')")
        }
        cli_alert("{.val {motif}}, score ({final_metric}): {.val {r$score}}")
        return(r)
    }, .parallel = parallel)

    scores <- sapply(res_s_list, function(x) x$score)

    res <- res_s_list[[which.max(scores)]]

    cli_alert_info("Best motif: {.val {res$seed_motif}}, score ({final_metric}): {.val {max(scores)}}")

    res$sample_idxs <- sample_idxs

    res$pred <- compute_pwm(sequences, res$pssm, res$spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect)
    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])]))
    }

    res$kmers <- cand_kmers

    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset, parallel = parallel)
    }

    cli_alert_success("Finished running regression. Consensus: {.val {res$consensus}}")

    if (is_binary_response(response)) {
        cli_alert_success("KS test D: {.val {round(res$ks$statistic, digits=4)}}, p-value: {.val {res$ks$p.value}}")
    } else {
        cli_alert_success("R^2: {.val {round(res$r2, digits=4)}}")
    }

    return(res)
}

get_cand_kmers <- function(sequences, response, kmer_length, min_gap, max_gap, min_kmer_cor, verbose, parallel = FALSE, max_cands = 10, ...) {
    all_kmers <- plyr::ldply(cli_progress_along(kmer_length), function(i) {
        screen_kmers(sequences, response, kmer_length = kmer_length[i], min_gap = 0, max_gap = max_gap, ...) %>%
            mutate(len = kmer_length[i], verbose = FALSE) %>%
            suppressMessages()
    }, .parallel = parallel)

    best_kmer <- all_kmers$kmer[which.max(abs(all_kmers$max_r2))] # return at least one kmer

    all_kmers <- all_kmers %>%
        # filter by correlation
        filter(sqrt(max_r2) > min_kmer_cor) %>%
        dplyr::distinct(kmer, .keep_all = TRUE)

    cands <- all_kmers %>%
        slice_max(n = min(nrow(all_kmers), max_cands), order_by = abs(max_r2)) %>%
        ungroup() %>%
        arrange(abs(max_r2)) %>%
        pull(kmer)

    dist_mat <- stringdist::stringdistmatrix(cands, cands, method = "osa", nthread = 1)
    for (i in 1:nrow(dist_mat)) {
        if (min(dist_mat[i, ], na.rm = TRUE) < 2) {
            cands <- cands[-i]
            dist_mat[i, ] <- NA
            dist_mat[, i] <- NA
        }
    }

    cands <- unique(c(best_kmer, cands))

    return(cands)
}
