
regress_multiple_motifs <- function(sequences,
                                    response,
                                    motif,
                                    motif_length,
                                    score_metric,
                                    bidirect,
                                    spat_min,
                                    spat_max,
                                    spat_bin,
                                    min_nuc_prob,
                                    is_train,
                                    include_response,
                                    seed,
                                    verbose,
                                    kmer_length,
                                    motif_num = 2,
                                    ...) {
    cli_alert_info("Running regression of {.val {motif_num}} motifs")

    res <- list()
    resp <- response

    for (i in 1:motif_num) {
        cli_alert_info("Running iteration {.val {i}} of {.val {max_iter}}")
        cli_alert_info("motif: {.val {motif}}")
        reg_result <- regress_pwm(
            sequences,
            resp,
            motif = motif,
            motif_length = motif_length,
            bidirect = bidirect,
            spat_min = spat_min,
            spat_max = spat_max,
            spat_bin = spat_bin,
            min_nuc_prob = min_nuc_prob,
            is_train = is_train,
            include_response = include_response,
            seed = seed,
            verbose = verbose,
            kmer_length = kmer_length,
            multi_step = FALSE,
            ...
        )

        if (i == 1) {
            reg_result$model <- lm(response ~ reg_result$pred)
        } else {
            pred_mat <- cbind(sapply(res, function(x) x$pred), reg_result$pred)
            reg_result$model <- lm(response ~ pred_mat)
        }
        cli_alert_info("Current r^2: {.val {summary(reg_result$model)$r.squared}}")

        res[[i]] <- reg_result
        resp <- reg_result$model$resid
        motif <- NULL
    }

    # final_res <- list(
    #     pssm = reg_result$pssm,
    #     spat = reg_result$spat,
    #     pred = predict(reg_result$model, newdata = sequences)
    # )
}
