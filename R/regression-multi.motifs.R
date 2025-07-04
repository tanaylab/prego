#' @export
#' @rdname regress_pwm
regress_multiple_motifs <- function(sequences,
                                    response,
                                    motif = NULL,
                                    motif_length = 15,
                                    score_metric = "r2",
                                    bidirect = TRUE,
                                    spat_bin_size = NULL,
                                    spat_num_bins = NULL,
                                    spat_model = NULL,
                                    improve_epsilon = 0.0001,
                                    min_nuc_prob = 0.001,
                                    unif_prior = 0.05,
                                    include_response = TRUE,
                                    seed = 60427,
                                    verbose = FALSE,
                                    kmer_length = 8,
                                    multi_kmers = FALSE,
                                    final_metric = NULL,
                                    max_cands = 10,
                                    min_gap = 0,
                                    max_gap = 1,
                                    min_kmer_cor = 0.08,
                                    motif_num = 2,
                                    smooth_k = 100,
                                    consensus_single_thresh = 0.6,
                                    consensus_double_thresh = 0.85,
                                    internal_num_folds = 1,
                                    match_with_db = TRUE,
                                    motif_dataset = all_motif_datasets(),
                                    parallel = getOption("prego.parallel", FALSE),
                                    alternative = "less",
                                    sample_for_kmers = FALSE,
                                    sample_frac = NULL,
                                    sample_idxs = NULL,
                                    sample_ratio = 1,
                                    log_energy = FALSE,
                                    energy_func_generator = NULL,
                                    energy_func = NULL,
                                    optimize_pwm = TRUE,
                                    optimize_spat = TRUE,
                                    ...) {
    if (motif_num < 2) {
        cli_abort("{.field motif_num} must be at least 2")
    }

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    max_seq_len <- nchar(sequences[1])
    bins <- calculate_bins(max_seq_len, spat_num_bins, spat_bin_size)
    spat_num_bins <- bins$spat_num_bins
    spat_bin_size <- bins$spat_bin_size
    cli_alert_info("Using {.val {spat_num_bins}} bins of size {.val {spat_bin_size}} bp")

    regression_func <- purrr::partial(regress_pwm,
        sequences = sequences,
        motif_length = motif_length,
        bidirect = bidirect,
        spat_bin_size = spat_bin_size,
        spat_num_bins = spat_num_bins,
        improve_epsilon = improve_epsilon,
        min_nuc_prob = min_nuc_prob,
        unif_prior = unif_prior,
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
        sample_for_kmers = sample_for_kmers,
        sample_frac = sample_frac,
        sample_idxs = sample_idxs,
        sample_ratio = sample_ratio,
        log_energy = log_energy,
        motif_num = 1,
        optimize_pwm = optimize_pwm,
        optimize_spat = optimize_spat,
        ...
    )

    cli_alert_info("Running regression for {.val {motif_num}} motifs")
    models <- list()
    e <- list()
    r0 <- response
    scores <- c()
    comb_scores <- c()

    cli_h2("Running first regression")
    res <- regression_func(response = response, motif = motif, score_metric = score_metric, final_metric = final_metric, spat_model = spat_model, energy_func = energy_func)

    if (!is.null(energy_func_generator)) {
        res_efunc <- apply_energy_func(
            prev_reg = res,
            response = response,
            regression_func = regression_func,
            score_metric = score_metric,
            final_metric = final_metric,
            energy_func_generator = energy_func_generator
        )
        cli::cli_alert_info("The second run changed the score from {.val {res$score}} to {.val {res_efunc$score}}")
        if (res_efunc$score > res$score) {
            res <- res_efunc
            cli::cli_alert_info("Taking the second run")
        } else {
            cli::cli_alert_info("Taking the first run")
        }
    }
    models[[1]] <- res
    e[[1]] <- res$pred
    names(e)[1] <- "e1"
    e_comb <- e[[1]]

    if (is_binary_response(response)) {
        comb_scores[1] <- res$ks$statistic
    } else {
        comb_scores[1] <- mean(res$r2)
    }
    scores[1] <- comb_scores[1]

    for (i in 2:motif_num) {
        cli_h2("Running regression #{i}")
        r <- r0 - pred_r_given_e(e_comb, r0, k = smooth_k) # residual
        res <- regression_func(response = r, score_metric = "r2", final_metric = "r2", spat_model = spat_model, energy_func = energy_func)
        if (!is.null(energy_func_generator)) {
            res_efunc <- apply_energy_func(
                prev_reg = res,
                response = response,
                regression_func = regression_func,
                score_metric = "r2",
                final_metric = "r2",
                energy_func_generator = energy_func_generator
            )
            cli::cli_alert_info("The second run changed the score from {.val {res$score}} to {.val {res_efunc$score}}")
            if (res_efunc$score > res$score) {
                res <- res_efunc
                cli::cli_alert_info("Taking the second run")
            } else {
                cli::cli_alert_info("Taking the first run")
            }
        }
        models[[i]] <- res

        e[[i]] <- res$pred
        names(e)[i] <- paste0("e", i)
        # combined lm model
        pred_df <- cbind(r0, as.data.frame(e))

        # Build a combined linear model. When the response has multiple columns we fit a
        # multivariate linear model using cbind(resp1, resp2, ...) ~ predictors.  When the
        # response is a single vector we simply model it as response ~ predictors.
        if (is.null(dim(r0)) || ncol(as.matrix(r0)) == 1) {
            # Single response variable
            pred_df <- as.data.frame(pred_df)
            colnames(pred_df)[1] <- "response"
            model_comb <- lm(response ~ ., data = pred_df)
            e_comb <- predict(model_comb, pred_df)
        } else {
            # Multi-response case
            r_cols <- colnames(r0)
            if (is.null(r_cols)) {
                # Provide default column names if missing
                r_cols <- paste0("resp", seq_len(ncol(r0)))
                colnames(pred_df)[seq_len(ncol(r0))] <- r_cols
            }
            # Construct the multivariate formula cbind(col1, col2, ...) ~ .
            resp_form <- paste0("cbind(", paste(r_cols, collapse = ", "), ")")
            model_comb <- lm(as.formula(paste(resp_form, "~ .")), data = pred_df)
            e_comb <- predict(model_comb, pred_df)
        }

        if (is_binary_response(response)) {
            ks <- suppressWarnings(ks.test(e[[i]][r0 == 1], e[[i]][r0 == 0], alternative = alternative)$statistic)
            cli_alert_info("KS statistic: {.val {ks}}")
            ks_comb <- suppressWarnings(ks.test(e_comb[r0 == 1], e_comb[r0 == 0], alternative = alternative)$statistic)
            cli_alert_info("KS test statistic for models {.val {1:i}}: {.val {ks_comb}}")
            cli_alert_info("Improvement in KS test statistic: {.val {ks_comb - comb_scores[i - 1]}}")
            comb_scores <- c(comb_scores, ks_comb)
            scores <- c(scores, ks)
        } else {
            r2 <- mean(cor(e[[i]], r0)^2)
            if (is.matrix(r0)) {
                r2_comb <- mean(diag(cor(e_comb, r0))^2)
            } else {
                r2_comb <- mean(cor(e_comb, r0)^2)
            }
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
            db_match_cor = sapply(models, function(x) mean(x$db_match_cor)),
            db_match_r2 = sapply(models, function(x) mean(x$db_match_r2))
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

    res$predict_multi <- function(x, parallel = getOption("prego.parallel", FALSE)) {
        e <- safe_llply(models, function(.x) .x$predict(x), .parallel = parallel) %>%
            do.call(cbind, .) %>%
            as.data.frame()
        colnames(e) <- paste0("e", seq_along(models))
        return(e)
    }

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[response == 1], res$pred[response == 0], alternative = alternative))
    } else {
        res$r2 <- cor(res$pred, response)^2
    }

    spat <- calc_spat_min_max(spat_num_bins, nchar(sequences[1]), spat_bin_size)

    res$spat_min <- spat$spat_min
    res$spat_max <- spat$spat_max
    res$spat_bin_size <- spat_bin_size
    res$bidirect <- bidirect
    res$seq_length <- nchar(sequences[1])

    res$export <- function(fn) {
        export_multi_regression(res, fn)
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

apply_energy_func <- function(prev_reg, response, regression_func, score_metric, final_metric, energy_func_generator) {
    energy_func <- energy_func_generator(prev_reg, response)
    cli::cli_alert_info("Running a second regression with energy function, initialized by the first regression")
    regression_func(response = response, motif = prev_reg$pssm, spat_model = prev_reg$spat, score_metric = score_metric, final_metric = final_metric, energy_func = energy_func, energy_func_generator = NULL)
}
