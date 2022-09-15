#' Run PWM regression on clusters
#'
#' @param clusters a vector with the cluster assignments for each sequence
#' @param two_phase whether to use two-phase optimization or not.
#' @param two_phase_sample_frac a vector of two numbers, specifying the fraction of
#' sequences to use in the first phase of optimization for the sequences which are not
#' in the cluster (first number) and in the cluster (second number).
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{models: }{a list with the models for each cluster}
#' \item{cluster_mat: }{an indicator matrix with the cluster assignments}
#' \item{pred_mat: }{a matrix with the predicted pwm for each sequence (rows) and cluster (columns)}
#' \item{stats: }{a data frame with the statistics for each cluster}
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
regress_pwm.clusters <- function(sequences, clusters, two_phase = TRUE, two_phase_sample_frac = c(0.1, 1), first_phase_metric = "ks", parallel = getOption("prego.parallel", TRUE), ...) {
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

    cli_alert_info("Running regression for {.val {nrow(cluster_mat)}} clusters")
    cluster_models <- plyr::llply(seq_len(ncol(cluster_mat)), function(i) {
        cli_h1("Cluster {.val {i}}")
        regression_func(sequences, cluster_mat[, i], ...)
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

    res <- list(
        models = cluster_models,
        cluster_mat = cluster_mat,
        pred_mat = pred_mat,
        stats = stats
    )

    return(res)
}
