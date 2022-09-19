#' Run PWM regression on a sample of the data
#'
#' @description The optimization would be performed with a sampled dataset of size \code{sample_frac}, or explicit sampled indices \code{sample_idxs}. Note that \code{multi_kmers} is TRUE by default.
#'
#' @param sample_frac fraction of the dataset to sample. When \code{response} is categorical (0 and 1), the sampling would be stratified by the category, i.e. \code{sample_frac} can be a vector of length 2 with the fraction of 0 and 1 responses to sample respectively.
#' If NULL - the default would be 0.1 for continuous variables, and for binary variables - the number of 0 responses would be equal to \code{sample_ratio} times the number of 1 responses.
#' @param sample_idxs indices of the sequences to use. If NULL, the indices would be sampled using \code{sample_frac}.
#' @param sample_ratio ratio between the '1' category and the '0' category in the sampled dataset. Relevant only when \code{sample_frac} is NULL.
#' @param verbose verbosity of the optimization.
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
                               kmer_length = 6:8,
                               multi_kmers = TRUE,
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

    cli_alert_info("Performing sampled optimization")
    if (is.null(sample_idxs)) {
        sample_idxs <- sample_response(response, sample_frac, sample_ratio, seed)
    }

    sequences_s <- sequences[sample_idxs]
    response_s <- response[sample_idxs, , drop = FALSE]

    res <- regress_pwm(
        sequences = sequences_s,
        response = response_s,
        motif = motif,
        multi_kmers = multi_kmers,
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
        kmer_length = kmer_length,
        max_cands = max_cands,
        min_gap = min_gap,
        max_gap = max_gap,
        min_kmer_cor = min_kmer_cor,
        consensus_single_thresh = consensus_single_thresh,
        consensus_double_thresh = consensus_double_thresh,
        final_metric = final_metric,
        match_with_db = FALSE,
        parallel = parallel
    )

    res$sample_idxs <- sample_idxs

    res$pred <- compute_pwm(sequences, res$pssm, res$spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect)
    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])]))
    }

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
        arrange(desc(abs(max_r2)))


    dist_mat <- stringdist::stringdistmatrix(cands$kmer, cands$kmer, method = "osa", nthread = 1)
    dist_mat[dist_mat != 1] <- NA
    g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected")
    cands <- cands %>%
        mutate(kmer_clust = igraph::cluster_louvain(g)$membership) %>%
        group_by(kmer_clust) %>%
        slice(1) %>%
        pull(kmer)

    cands <- unique(c(best_kmer, cands))

    return(cands)
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
        sample_idxs <- sample(seq_along(sequences), size = floor(sample_frac * length(sequences)))
    }

    if (is_binary_response(response)) {
        cli_alert_info("sampled {.val {sum(response[sample_idxs, 1] == 0)}} zeros and {.val {sum(response[sample_idxs, 1] == 1)}} ones")
    }
    return(sample_idxs)
}
