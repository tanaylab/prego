
#' @export
#' @rdname regress_pwm
regress_multiple_motifs <- function(sequences,
                                    response,
                                    motif = NULL,
                                    motif_length = 15,
                                    score_metric = "r2",
                                    bidirect = TRUE,
                                    spat_min = 0,
                                    spat_max = NULL,
                                    spat_bin = 50,
                                    spat_model = NULL,
                                    improve_epsilon = 0.0001,
                                    min_nuc_prob = 0.001,
                                    unif_prior = 0.05,
                                    is_train = NULL,
                                    include_response = TRUE,
                                    seed = 60427,
                                    verbose = FALSE,
                                    kmer_length = 8,
                                    multi_kmers = FALSE,
                                    final_metric = NULL,
                                    max_cands = 10,
                                    min_gap = 0,
                                    max_gap = 1,
                                    min_kmer_cor = 0.1,
                                    motif_num = 2,
                                    smooth_k = 100,
                                    consensus_single_thresh = 0.6,
                                    consensus_double_thresh = 0.85,
                                    internal_num_folds = 1,
                                    match_with_db = TRUE,
                                    motif_dataset = all_motif_datasets(),
                                    parallel = getOption("prego.parallel", FALSE),
                                    alternative = "less",
                                    ...) {
    if (motif_num < 2) {
        cli_abort("{.field motif_num} must be at least 2")
    }

    regression_func <- purrr::partial(regress_pwm,
        sequences = sequences,
        motif_length = motif_length,
        bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_bin = spat_bin,
        spat_model = spat_model,
        improve_epsilon = improve_epsilon,
        min_nuc_prob = min_nuc_prob,
        unif_prior = unif_prior,
        is_train = is_train,
        multi_kmers = multi_kmers,
        include_response = include_response,
        seed = seed,
        verbose = verbose,
        kmer_length = kmer_length,
        max_cands = max_cands,
        min_gap = min_gap,
        max_gap = max_gap,
        min_kmer_cor = min_kmer_cor,
        consensus_single_thresh = consensus_single_thresh,
        consensus_double_thresh = consensus_double_thresh,
        internal_num_folds = internal_num_folds,
        parallel = parallel,
        match_with_db = match_with_db,
        motif_dataset = motif_dataset,
        alternative = alternative,
        ...
    )

    cli_alert_info("Running regression for {.val {motif_num}} motifs")
    models <- list()
    e <- list()
    r0 <- response
    scores <- c()
    comb_scores <- c()

    cli_h2("Running first regression")
    res <- regression_func(response = response, motif = motif, score_metric = score_metric, final_metric = final_metric)
    models[[1]] <- res
    e[[1]] <- res$pred
    names(e)[1] <- "e1"
    e_comb <- e[[1]]

    if (is_binary_response(response)) {
        comb_scores[1] <- res$ks$statistic
    } else {
        comb_scores[1] <- res$r2
    }
    scores[1] <- comb_scores[1]

    for (i in 2:motif_num) {
        cli_h2("Running regression #{i}")
        r <- r0 - pred_r_given_e(e_comb, r0, k = smooth_k) # residual
        res <- regression_func(response = r, score_metric = "r2", final_metric = "r2")
        models[[i]] <- res

        e[[i]] <- res$pred
        names(e)[i] <- paste0("e", i)

        # combined lm model
        pred_df <- cbind(r0, as.data.frame(e))
        model_comb <- lm(r0 ~ ., data = pred_df)
        e_comb <- predict(model_comb, pred_df)

        if (is_binary_response(response)) {
            ks <- suppressWarnings(ks.test(e[[i]][r0 == 1], e[[i]][r0 == 0], alternative = alternative)$statistic)
            cli_alert_info("KS statistic: {.val {ks}}")
            ks_comb <- suppressWarnings(ks.test(e_comb[r0 == 1], e_comb[r0 == 0], alternative = alternative)$statistic)
            cli_alert_info("KS test statistic for models {.val {1:i}}: {.val {ks_comb}}")
            cli_alert_info("Improvement in KS test statistic: {.val {ks_comb - comb_scores[i - 1]}}")
            comb_scores <- c(comb_scores, ks_comb)
            scores <- c(scores, ks)
        } else {
            r2 <- cor(e[[i]], r0)^2
            r2_comb <- cor(e_comb, r0)^2
            cli_alert_info("R2 for models {.val {1:i}}: {.val {r2_comb}}")
            cli_alert_info("Improvement in R2: {.val {r2_comb - comb_scores[i - 1]}}")
            comb_scores <- c(comb_scores, r2_comb)
            scores <- c(scores, r2)
        }
    }

    stats <- data.frame(
        model = 1:motif_num,
        score = scores,
        comb_score = comb_scores,
        diff = c(NA, diff(comb_scores)),
        consensus = sapply(models, function(x) x$consensus),
        seed_motif = sapply(models, function(x) x$seed_motif)
    )

    if (match_with_db) {
        stats <- stats %>% mutate(
            db_match = sapply(models, function(x) x$db_match),
            db_match_cor = sapply(models, function(x) x$db_match_cor),
            db_match_r2 = sapply(models, function(x) x$db_match_r2)
        )
        if (is_binary_response(response)) {
            stats <- stats %>% mutate(
                db_match_ks = sapply(models, function(x) x$db_match_ks$statistic)
            )
        }
    }

    best_model <- which.max(scores)
    cli_alert_success("Best model: model #{best_model} (score of {.val {scores[best_model]}})")
    res <- list(
        models = models,
        multi_stats = stats,
        model = model_comb # last combined model
    )
    res$predict <- function(x, ...) {
        e <- lapply(1:motif_num, function(i) models[[i]]$predict(x, ...))
        e <- as.data.frame(e)
        colnames(e) <- paste0("e", 1:motif_num)
        predict(model_comb, e, ...)
    }
    res$pred <- res$predict(sequences)

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[response == 1], res$pred[response == 0], alternative = alternative))
    } else {
        res$r2 <- cor(res$pred, response)^2
    }

    return(res)
}

pred_r_given_e <- function(e, r, k = 100) {
    k <- min(k, length(r))
    tibble(e = e, r = r) %>%
        mutate(i = 1:n()) %>%
        arrange(e) %>%
        mutate(pred = zoo::rollapply(r, k, mean, partial = TRUE)) %>%
        arrange(i) %>%
        pull(pred)
}
