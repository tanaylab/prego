
regress_multiple_motifs <- function(sequences,
                                    response,
                                    motif,
                                    motif_length,
                                    score_metric,
                                    bidirect,
                                    spat_min,
                                    spat_max,
                                    spat_bin,
                                    spat_model,
                                    improve_epsilon,
                                    min_nuc_prob,
                                    unif_prior,
                                    is_train,
                                    include_response,
                                    seed,
                                    verbose,
                                    kmer_length,
                                    motif_num,
                                    consensus_single_thresh,
                                    consensus_double_thresh,
                                    match_with_db,
                                    ...) {
    cli_abort("Multiple motifs are not supported yet")
    # cli_alert_info("Running regression of {.val {motif_num}} motifs")
    # first_regression_func <- purrr::partial(
    #     regress_pwm,
    #     motif = motif,
    #     motif_length = motif_length,
    #     score_metric = score_metric,
    #     bidirect = bidirect,
    #     spat_min = spat_min,
    #     spat_max = spat_max,
    #     spat_bin = spat_bin,
    #     spat_model = spat_model,
    #     improve_epsilon = improve_epsilon,
    #     min_nuc_prob = min_nuc_prob,
    #     unif_prior = unif_prior,
    #     is_train = is_train,
    #     include_response = FALSE,
    #     seed = seed,
    #     verbose = verbose,
    #     kmer_length = kmer_length,
    #     consensus_single_thresh = consensus_single_thresh,
    #     consensus_double_thresh = consensus_double_thresh,
    #     match_with_db = match_with_db,
    #     ...
    # )


    # res <- list()
    # resp <- response

    # for (i in 1:motif_num) {
    #     cli_alert_info("Running iteration {.val {i}} of {.val {max_iter}}")
    #     cli_alert_info("motif: {.val {motif}}")
    #     reg_result <- regress_pwm(
    #             sequences = sequences,
    #             response = response,
    #             motif = NULL,
    #             motif_length = motif_length,
    #             score_metric = score_metric,
    #             bidirect = bidirect,
    #             spat_min = spat_min,
    #             spat_max = spat_max,
    #             spat_bin = spat_bin,
    #             improve_epsilon = improve_epsilon,
    #             min_nuc_prob = min_nuc_prob,
    #             unif_prior = unif_prior,
    #             is_train = is_train,
    #             include_response = include_response,
    #             seed = seed,
    #             verbose = verbose,
    #             kmer_length = kmer_length,
    #             motif_num = 1,
    #             ...
    #     )

    #     if (i == 1) {
    #         reg_result$model <- lm(response ~ reg_result$pred)
    #     } else {
    #         pred_mat <- cbind(sapply(res, function(x) x$pred), reg_result$pred)
    #         reg_result$model <- lm(response ~ pred_mat)
    #     }
    #     cli_alert_info("Current r^2: {.val {summary(reg_result$model)$r.squared}}")

    #     res[[i]] <- reg_result
    #     resp <- reg_result$model$resid
    #     motif <- NULL
    # }

    # # final_res <- list(
    # #     pssm = reg_result$pssm,
    # #     spat = reg_result$spat,
    # #     pred = predict(reg_result$model, newdata = sequences)
    # # )
}
