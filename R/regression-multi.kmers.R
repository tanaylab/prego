regress_pwm.multi_kmers <- function(sequences,
                                    response,
                                    motif_length = 15,
                                    score_metric = "r2",
                                    bidirect = TRUE,
                                    spat_bin_size = NULL,
                                    spat_num_bins = 7,
                                    spat_model = NULL,
                                    improve_epsilon = 0.0001,
                                    min_nuc_prob = 0.001,
                                    unif_prior = 0.05,
                                    include_response = TRUE,
                                    seed = 60427,
                                    verbose = FALSE,
                                    kmer_length = 6:8,
                                    max_cands = 10,
                                    min_gap = 0,
                                    max_gap = 1,
                                    min_kmer_cor = 0.08,
                                    val_frac = 0.1,
                                    consensus_single_thresh = 0.6,
                                    consensus_double_thresh = 0.85,
                                    internal_num_folds = 1,
                                    final_metric = "r2",
                                    parallel = getOption("prego.parallel", FALSE),
                                    match_with_db = TRUE,
                                    screen_db = FALSE,
                                    motif_dataset = all_motif_datasets(),
                                    alternative = "two.sided",
                                    sample_for_kmers = FALSE,
                                    sample_frac = NULL,
                                    sample_idxs = NULL,
                                    sample_ratio = 1,
                                    log_energy = FALSE,
                                    energy_func = NULL,
                                    xmin = -100,
                                    xmax = 100,
                                    npts = 1e4,
                                    optimize_pwm = TRUE,
                                    optimize_spat = TRUE,
                                    kmer_sequence_length = NULL,
                                    symmetrize_spat = TRUE,
                                    ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(dim(response))) {
        response <- matrix(response, ncol = 1)
    } else if (!is.matrix(response)) {
        response <- as.matrix(response)
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    if (sample_for_kmers) {
        cli_alert_info("Performing sampled optimization")
        if (is.null(sample_idxs)) {
            sample_idxs <- sample_response(response, sample_frac, sample_ratio, seed)
        }

        sequences_s <- sequences[sample_idxs]
        response_s <- response[sample_idxs, , drop = FALSE]
    } else {
        sequences_s <- sequences
        response_s <- response
    }

    # split to validation and train
    val_idxs <- sample_response(response_s, val_frac, 1, seed)
    train_idxs <- setdiff(1:length(sequences_s), val_idxs)

    regress_pwm_single_kmer <- purrr::partial(
        regress_pwm,
        motif_length = motif_length,
        score_metric = score_metric,
        bidirect = bidirect,
        spat_bin_size = spat_bin_size,
        spat_num_bins = spat_num_bins,
        improve_epsilon = improve_epsilon,
        min_nuc_prob = min_nuc_prob,
        unif_prior = unif_prior,
        include_response = FALSE,
        seed = seed,
        verbose = FALSE,
        consensus_single_thresh = consensus_single_thresh,
        consensus_double_thresh = consensus_double_thresh,
        internal_num_folds = internal_num_folds,
        match_with_db = FALSE,
        alternative = alternative,
        multi_kmers = FALSE,
        sample_for_kmers = FALSE,
        log_energy = log_energy,
        energy_func = energy_func,
        xmin = xmin,
        xmax = xmax,
        npts = npts,
        optimize_pwm = optimize_pwm,
        optimize_spat = optimize_spat
    )

    if (!is.null(kmer_sequence_length)) {
        if (kmer_sequence_length > nchar(sequences[1])) {
            cli_abort("kmer_sequence_length cannot be greater than the length of the sequences")
        }
        # str_sub the sequence from the middle
        sequences_kmers <- stringr::str_sub(sequences, start = (nchar(sequences) - kmer_sequence_length) / 2 + 1, end = (nchar(sequences) + kmer_sequence_length) / 2)
        cli_alert_info("Using kmer_sequence_length of {.val {kmer_sequence_length}}")
    } else {
        sequences_kmers <- sequences
    }

    cli_h3("Generate candidate kmers")
    cand_kmers <- get_cand_kmers(sequences_kmers, response, kmer_length, min_gap, max_gap, min_kmer_cor, verbose, parallel, max_cands = max_cands, ...)

    if (sample_for_kmers) {
        cli_h3("Regress each candidate kmer on sampled data")
    } else {
        cli_h3("Regress each candidate kmer")
    }

    cli_alert_info("Running regression on {.val {length(cand_kmers)}} candidate kmers")
    cli_ul(c(
        "Bidirectional: {.val {bidirect}}",
        "Spat bin size: {.val {spat_bin_size}}",
        "Number of spatial bins {.val {spat_num_bins}}",
        "Length of sequence: {.val {spat_bin_size*spat_num_bins}}",
        "Min gap: {.val {min_gap}}",
        "Max gap: {.val {max_gap}}",
        "Kmer length: {.val {kmer_length}}",
        "Improve epsilon: {.val {improve_epsilon}}",
        "Min nuc prob: {.val {min_nuc_prob}}",
        "Uniform prior: {.val {unif_prior}}",
        "Score metric: {.val {score_metric}}",
        "Seed: {.val {seed}}"
    ))

    run_kmer <- function(i) {
        motif <- cand_kmers[i]
        cli_alert("regressing with seed: {.val {motif}}")
        r <- regress_pwm_single_kmer(motif = motif, sequences = sequences_s[train_idxs], response = response_s[train_idxs]) %>%
            suppressMessages()
        pr <- r$predict(sequences_s[val_idxs])
        if (final_metric == "ks") {
            if (!is_binary_response(response_s)) {
                cli_abort("Cannot use {.field final_metric} {.val ks} when {.field response} is not binary")
            }
            r$score <- r[[final_metric]]$statistic
            r$val_score <- suppressWarnings(ks.test(pr[as.logical(response_s[val_idxs, 1])], pr[!as.logical(response_s[val_idxs, 1])], alternative = alternative)$statistic)
        } else if (final_metric == "r2") {
            r$score <- r[[final_metric]]
            r$val_score <- cor(pr, response_s[val_idxs])^2
        } else {
            cli_abort("Unknown {.field final_metric} (can be 'ks' or 'r2')")
        }
        cli_alert("{.val {motif}}, score ({final_metric}): {.val {r$score}}, validation score ({final_metric}): {.val {r$val_score}}")
        return(r)
    }

    res_kmer_list <- safe_llply(cli_progress_along(cand_kmers), run_kmer, .parallel = parallel)

    # find jobs that failed (not numeric)
    failed <- sapply(res_kmer_list, function(x) !is.numeric(x$score))

    if (any(failed)) {
        cli_alert_warning("Some kmers failed to regress")
        cli_alert_info("Performing regression on full data")
        res_kmer_list[failed] <- lapply(cand_kmers[failed], function(x) regress_pwm_single_kmer(motif = x, sequences = sequences, response = response))
    }

    scores <- sapply(res_kmer_list, function(x) x$val_score)
    if (is.matrix(scores) && nrow(scores) > 1) {
        scores <- colMeans(scores)
    }

    if (length(which.max(scores)) == 0) {
        cli_alert_warning("No motifs found")
        cli_alert_info("Performing regression on full data")
        res <- regress_pwm_single_kmer(motif = cand_kmers[1], sequences = sequences, response = response)
    } else {
        cli_alert_info("Performing regression on full data")
        res <- regress_pwm_single_kmer(motif = cand_kmers[which.max(scores)], sequences = sequences, response = response)
    }

    res$kmers <- cand_kmers

    if (include_response) {
        res$response <- response
    }

    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset, alternative = alternative, parallel = parallel)
    }

    if (screen_db) {
        res <- add_regression_db_screen(res, response, sequences, motif_dataset, final_metric, alternative = alternative, prior = unif_prior, bidirect = bidirect, parallel = parallel)
    }

    cli_alert_info("Best motif: {.val {res$seed_motif}}, score ({final_metric}): {.val {max(scores)}}")

    spat <- calc_spat_min_max(spat_num_bins, nchar(sequences[1]), spat_bin_size)

    res$spat_min <- spat$spat_min
    res$spat_max <- spat$spat_max
    res$spat_bin_size <- spat_bin_size
    res$bidirect <- bidirect
    res$seq_length <- nchar(sequences[1])

    res <- add_predict_function(res, spat, bidirect, energy_func)
    res$pred <- res$predict(sequences)

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2
    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])], alternative = alternative))
        res$score <- res$ks$statistic
    } else {
        res$score <- res$r2
    }

    return(res)
}

get_cand_kmers <- function(sequences, response, kmer_length, min_gap, max_gap, min_kmer_cor, verbose, parallel = FALSE, max_cands = 10, ...) {
    all_kmers <- safe_ldply(cli_progress_along(kmer_length), function(i) {
        screen_kmers(sequences, response, kmer_length = kmer_length[i], min_gap = min_gap, max_gap = max_gap, min_cor = min_kmer_cor, ...) %>%
            mutate(len = kmer_length[i], verbose = FALSE) %>%
            suppressMessages()
    }, .parallel = parallel)

    best_kmer <- all_kmers$kmer[which.max(abs(all_kmers$max_r2))] # return at least one kmer

    if (length(best_kmer) == 0) { # could not find any kmer
        cli::cli_alert_info("Could not find any kmer when using a threshold of {.val {min_kmer_cor}}. Trying with {.val {min_kmer_cor/2}}")
        all_kmers <- safe_ldply(cli_progress_along(kmer_length), function(i) {
            screen_kmers(sequences, response, kmer_length = kmer_length[i], min_gap = min_gap, max_gap = max_gap, min_cor = min_kmer_cor / 2, ...) %>%
                mutate(len = kmer_length[i], verbose = FALSE) %>%
                suppressMessages()
        }, .parallel = parallel)
        best_kmer <- all_kmers$kmer[which.max(abs(all_kmers$max_r2))]
        if (length(best_kmer) == 0) { # could not find any kmer
            res <- paste(rep("*", kmer_length), collapse = "")
            cli_alert_info("Could not find any kmer. Initializing with {.val {res}}")
            return(res)
        }
    }

    all_kmers <- all_kmers %>%
        # filter by correlation
        filter(sqrt(max_r2) > min_kmer_cor) %>%
        dplyr::distinct(kmer, .keep_all = TRUE)

    cands <- all_kmers %>%
        slice_max(n = min(nrow(all_kmers), max_cands), order_by = abs(max_r2)) %>%
        arrange(desc(abs(max_r2)))


    dist_mat <- stringdist::stringdistmatrix(cands$kmer, cands$kmer, method = "osa", nthread = 1)
    if (any(dist_mat == 1) && ncol(dist_mat) == nrow(dist_mat)) {
        dist_mat[dist_mat != 1] <- NA
        dist_mat[is.na(dist_mat)] <- 0
        g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected")
        cands <- cands %>%
            mutate(kmer_clust = igraph::cluster_louvain(g)$membership) %>%
            group_by(kmer_clust) %>%
            slice(1) %>%
            pull(kmer)
    } else {
        cands <- cands$kmer
    }

    cands <- unique(c(best_kmer, cands))

    return(cands)
}
