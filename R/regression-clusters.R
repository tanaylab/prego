#' Run PWM regression on clusters
#'
#' @param clusters a vector with the cluster assignments for each sequence
#' @param use_sample whether to use sampled optimization or not (default: FALSE).
#' @param sample_frac a vector of two numbers, specifying the fraction of
#' sequences to use in when sampling for the sequences which are not
#' in the cluster (first number) and in the cluster (second number). If NULL -
#' @param sample_ratio When \code{sample_frac} is NULL, the number of sequences not in the cluster would be equal to \code{sample_ratio} times the number of sequences in the cluster.
#' @param match_with_db match the resulting PWMs with motif databases using \code{pssm_match}.
#' This would add a column named 'db_match' to the stats data frame, together with 'pred_mat_db' with the
#' database motif predictions, and and 'db_dataset' which is similiar to 'motif_dataset' for the database motifs.
#' Note that the closest match is returned, even if it is not similar enough in absolute terms.
#' Also, the match is done between the rsulting regression \emph{pssm} and the pssms in the databse - in order to find the best motif in the database which explain the clusters, use \code{screen_pwm.clusters}.
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{models: }{a list with the models for each cluster}
#' \item{cluster_mat: }{an indicator matrix with the cluster assignments}
#' \item{pred_mat: }{a matrix with the predicted pwm for each sequence (rows) and cluster (columns)}
#' \item{motif_dataset: }{a data frame with the PSSMs for each cluster}
#' \item{spat_dataset: }{a data frame with the spatial model for each cluster}
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
#' # multiple motifs per cluster
#' res_multi <- regress_pwm.clusters(cluster_sequences_example, clusters_example, motif_num = 3)
#' res_multi$multi_stats
#' plot_regression_qc_multi(res_multi$models[[1]], title = names(res$models)[1])
#'
#' @inheritParams regress_pwm
#' @inheritParams regress_pwm.sample
#' @inheritDotParams regress_pwm
#' @inheritDotParams regress_pwm.sample
#'
#' @export
regress_pwm.clusters <- function(sequences, clusters, use_sample = TRUE, match_with_db = TRUE, sample_frac = NULL, sample_ratio = 1, final_metric = "ks", parallel = getOption("prego.parallel", TRUE), ...) {
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
    if (use_sample) {
        cli_alert_info("Using sampled optimization")
        regression_func <- purrr::partial(regress_pwm.sample, sample_frac = sample_frac, final_metric = final_metric, sample_ratio = sample_ratio, parallel = FALSE)
        if ("sample_idxs" %in% names(list(...))) {
            cli_abort("The {.field sample_idxs} argument is not supported in {.fun regress_pwm.clusters}")
        }
    } else {
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

    motif_dataset <- purrr::imap_dfr(cluster_models, ~ .x$pssm %>% mutate(motif = .y)) %>%
        select(motif, pos, A, C, G, T)
    spat_dataset <- purrr::imap_dfr(cluster_models, ~ .x$spat %>% mutate(motif = .y)) %>%
        select(motif, bin, spat_factor)

    res <- list(
        models = cluster_models,
        cluster_mat = cluster_mat,
        pred_mat = pred_mat,
        motif_dataset = motif_dataset,
        spat_dataset = spat_dataset,
        stats = stats
    )

    if (match_with_db) {
        cli_alert_info("Matching with motif databases")
        res$stats$db_match <- purrr::map_chr(cluster_models, "db_match")
        res$stats$db_match_dist <- purrr::map_chr(cluster_models, "db_match_dist")
        res$pred_mat_db <- purrr::map(cluster_models, "db_match_pred") %>% do.call(cbind, .)
        colnames(res$pred_mat_db) <- names(res$stats$db_match)
        rownames(res$pred_mat_db) <- names(sequences)
        res$db_dataset <- purrr::imap_dfr(cluster_models, ~ .x$db_match_pssm %>% mutate(cluster = .y, motif = .x$db_match)) %>%
            select(cluster, motif, pos, A, C, G, T)
    }

    if ("models" %in% names(res$models[[1]])) {
        res$multi_stats <- purrr::imap_dfr(cluster_models, ~ .x$multi_stats %>% mutate(cluster = .y)) %>%
            select(cluster, everything())
    }

    return(res)
}

#' Screen for motifs in a database for every cluster
#'
#' @param sequences a vector with the sequences
#' @param clusters a vector with the cluster assignments
#' @param min_D minimum distance to consider a match
#' @param only_match if TRUE, only return the best match for each cluster
#'
#' @return a matrix with the KS D statistics for each cluster (columns) and every motif (rows)
#' that had at least one cluster with D >= min_D. If \code{only_match} is TRUE, a named vector
#' with the name of best motif match for each cluster is returned (regardless of \code{min_D}).
#'
#' @examples
#' D_mat <- screen_pwm.clusters(cluster_sequences_example, clusters_example)
#' dim(D_mat)
#' D_mat[1:5, 1:5]
#'
#' # return only the best match
#' screen_pwm.clusters(cluster_sequences_example, clusters_example, only_match = TRUE)
#'
#' @inheritParams extract_pwm
#' @inheritDotParams compute_pwm
#' @export
screen_pwm.clusters <- function(sequences, clusters, dataset = all_motif_datasets(), motifs = NULL, parallel = getOption("prego.parallel", TRUE), min_D = 0.4, only_match = FALSE, prior = 0.01, ...) {
    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    sequences <- toupper(sequences)
    cluster_ids <- unique(clusters)

    res <- plyr::daply(dataset, "motif", function(x) {
        pwm <- compute_pwm(sequences, x, prior = prior, ...)
        purrr::map_dbl(cluster_ids, function(cl) {
            suppressWarnings(ks.test(pwm[clusters == cl], pwm[clusters != cl])$statistic)
        })
    }, .parallel = parallel)
    colnames(res) <- cluster_ids
    res <- as.matrix(res)

    if (only_match) {
        best_match <- rownames(res)[apply(res, 2, which.max)]
        names(best_match) <- cluster_ids
        return(best_match)
    }

    maxs <- apply(res, 1, max)
    f <- which(maxs >= min_D)
    if (length(f) == 0) {
        cli_abort("No matches found, try to lower the {.field min_D} parameter")
    }
    res <- res[f, , drop = FALSE]

    return(res)
}
