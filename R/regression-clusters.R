#' Run PWM regression on clusters
#'
#' @param clusters a vector with the cluster assignments for each sequence
#' @param two_phase whether to use two-phase optimization or not (default: FALSE).
#' @param two_phase_sample_frac a vector of two numbers, specifying the fraction of
#' sequences to use in the first phase of optimization for the sequences which are not
#' in the cluster (first number) and in the cluster (second number).
#' @param match_with_db match the resulting PWMs with motif databases using \code{pssm_match}.
#' This would add a column named 'db_match' to the stats data frame. Note that the closest match
#' is returned, even if it is not similar enough in absolute terms.
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{models: }{a list with the models for each cluster}
#' \item{cluster_mat: }{an indicator matrix with the cluster assignments}
#' \item{pred_mat: }{a matrix with the predicted pwm for each sequence (rows) and cluster (columns)}
#' \item{stats: }{a data frame with statistics for each cluster}
#' }
#'
#' @examples
#' res <- regress_pwm.clusters(cluster_sequences_example, clusters_example)
#' head(res$pred_mat)
#' res$stats
#'
#' plot_regression_qc(res$models[[1]], title = names(res$models)[1])
#'
#' @inheritParams regress_pwm
#' @inheritParams regress_pwm.two_phase
#' @inheritDotParams regress_pwm
#' @inheritDotParams regress_pwm.two_phase
#'
#' @export
regress_pwm.clusters <- function(sequences, clusters, two_phase = FALSE, match_with_db = TRUE, two_phase_sample_frac = c(0.1, 1), first_phase_metric = "ks", parallel = getOption("prego.parallel", TRUE), ...) {
    if (length(clusters) != length(sequences)) {
        cli_abort("The {.field clusters} vector should have the same length as the {.field sequences} vector")
    }

    na_idx <- which(is.na(clusters))
    if (length(na_idx) > 0) {
        cli_alert_warning("Removing {.val {length(na_idx)}} sequences with NA cluster assignments")
        sequences <- sequences[-na_idx]
        clusters <- clusters[-na_idx]
    }

    if (is.null(names(clusters))) {
        names(clusters) <- seq_along(clusters)
    }
    names(sequences) <- names(clusters)

    cluster_mat <- tibble::enframe(clusters, "id", "cluster") %>%
        mutate(i = 1) %>%
        tidyr::spread("cluster", i, fill = 0) %>%
        tibble::column_to_rownames("id") %>%
        as.matrix()

    cluster_mat <- cluster_mat[names(sequences), ]

    if (two_phase) {
        cli_alert_info("Using two-phase optimization")
        regression_func <- purrr::partial(regress_pwm.two_phase, two_phase_sample_frac = two_phase_sample_frac, first_phase_metric = first_phase_metric, parallel = FALSE)
        if ("first_phase_idxs" %in% names(list(...))) {
            cli_abort("The {.field first_phase_idxs} argument is not supported in {.fun regress_pwm.clusters}")
        }
    } else {
        cli_alert_info("Using single-phase optimization")
        regression_func <- regress_pwm
    }

    cli_alert_info("Running regression for {.val {ncol(cluster_mat)}} clusters")
    cluster_models <- plyr::llply(seq_len(ncol(cluster_mat)), function(i) {
        cli_h1("Cluster {.val {i}}")
        regression_func(sequences, cluster_mat[, i], match_with_db = match_with_db, ...)
    }, .parallel = parallel)

    names(cluster_models) <- colnames(cluster_mat)

    pred_mat <- purrr::map(cluster_models, "pred") %>% do.call(cbind, .)
    colnames(pred_mat) <- names(cluster_models)
    rownames(pred_mat) <- names(sequences)

    stopifnot(all.equal(names(cluster_models), colnames(pred_mat)))

    stats <- purrr::imap_dfr(cluster_models, ~ {
        tibble(
            cluster = .y,
            consensus = .x$consensus,
            ks_D = .x$ks$statistic,
            r2 = .x$r2,
            seed_motif = .x$seed_motif
        )
    })

    if (match_with_db) {
        cli_alert_info("Matching with motif databases")
        stats$db_match <- purrr::map_chr(cluster_models, "db_match")
        stats$db_match_dist <- purrr::map_chr(cluster_models, "db_match_dist")
        pred_mat_db <- purrr::map(cluster_models, "db_match_pred") %>% do.call(cbind, .)
        colnames(pred_mat_db) <- names(stats$db_match)
        rownames(pred_mat_db) <- names(sequences)
    }

    res <- list(
        models = cluster_models,
        cluster_mat = cluster_mat,
        pred_mat = pred_mat,
        stats = stats
    )

    if (match_with_db) {
        res$db_pred_mat <- pred_mat_db
    }

    return(res)
}
