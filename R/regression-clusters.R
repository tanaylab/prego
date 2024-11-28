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
#' database motif predictions, and 'db_dataset' which is similar to 'motif_dataset' for the database motifs.
#' Note that the closest match is returned, even if it is not similar enough in absolute terms.
#' Also, the match is done between the resulting regression \emph{pssm} and the pssms in the database - in order to find the best motif in the database set \code{screen_db=TRUE}.
#' @param screen_db screen for the best motif in the database which explains the clusters. See \code{screen_pwm.clusters}.
#' @param use_sge use the function \code{gcluster.run2} from the misha.ext package to run the optimization on a SGE cluster. Only relevant if the \code{misha.ext} package is installed. Note that \code{gcluster.run2} writes the current
#' environment before starting the parallelization, so it is better to run this function in a clean environment.
#' Also, Note that 'prego' needs to be installed in order for this to work, i.e. you cannot use \code{devtools::load_all()} or \code{pkgload::load_all()} to load the package.
#' @param alternative alternative hypothesis for the KS test. Can be "two.sided", "less" or "greater"
#'
#' @return a list with the following elements:
#' \describe{
#' \item{models: }{a list with the models for each cluster}
#' \item{cluster_mat: }{an indicator matrix with the cluster assignments}
#' \item{pred_mat: }{a matrix of the energies of the predicted motifs per cluster (columns) in each sequence (rows)}
#' \item{motif_dataset: }{a data frame with the PSSMs for each cluster}
#' \item{spat_dataset: }{a data frame with the spatial model for each cluster}
#' \item{stats: }{a data frame with statistics for each cluster}
#' }
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm.clusters(cluster_sequences_example, clusters_example)
#' head(res$pred_mat)
#' res$stats
#'
#' plot_regression_qc(res$models[[1]], title = names(res$models)[1])
#'
#' # multiple motifs per cluster
#' res_multi <- regress_pwm.clusters(cluster_sequences_example, clusters_example, motif_num = 3)
#' res_multi$multi_stats
#' plot_regression_qc_multi(res_multi$models[[1]], title = names(res_multi$models)[1])
#' }
#'
#' # screen also for the best motif in the database
#' res_screen <- regress_pwm.clusters(cluster_sequences_example, clusters_example, screen_db = TRUE)
#' res_screen$stats
#'
#' plot_regression_qc(res_screen$models[[1]], title = names(res_screen$models)[1])
#'
#' @inheritParams regress_pwm
#' @inheritParams regress_pwm.sample
#' @inheritDotParams regress_pwm
#' @inheritDotParams regress_pwm.sample
#' @inheritParams screen_pwm.clusters
#'
#' @export
regress_pwm.clusters <- function(sequences, clusters, use_sample = TRUE, match_with_db = TRUE, screen_db = FALSE, sample_frac = NULL, sample_ratio = 1, final_metric = "ks", parallel = getOption("prego.parallel", TRUE), use_sge = FALSE, dataset = all_motif_datasets(), motifs = NULL, min_D = 0, prior = 0.01, alternative = "two.sided", ...) {
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
        regression_func <- purrr::partial(regress_pwm.sample, sample_frac = sample_frac, final_metric = final_metric, sample_ratio = sample_ratio)
        if ("sample_idxs" %in% names(list(...))) {
            cli_abort("The {.field sample_idxs} argument is not supported in {.fun regress_pwm.clusters}")
        }
    } else {
        regression_func <- regress_pwm
    }

    cli_alert_info("Running regression for {.val {ncol(cluster_mat)}} clusters")
    if (use_sge) {
        if (!("misha.ext" %in% utils::installed.packages())) {
            cli_abort("The {.field misha.ext} package is required when {.code use_sge=TRUE}. Please install it with {.code remotes::install_packages('tanaylab/misha.ext')}.")
        }
        if (!("prego" %in% utils::installed.packages())) {
            cli_abort("The {.field prego} package needs to be installed when {.code use_sge=TRUE}. Please install it with {.code remotes::install_packages('tanaylab/prego')}.")
        }
        cli_alert_info("Using SGE cluster")
        cmds <- paste0("regression_func(sequences, cluster_mat[, ", seq_len(ncol(cluster_mat)), "], match_with_db = match_with_db, parallel = parallel, ...)")
        sge_res <- misha.ext::gcluster.run2(command_list = cmds)
        ret_class <- purrr::map_chr(sge_res, ~ class(.x$retv))
        if (any(ret_class != "list")) {
            failed <- which(ret_class != "list")
            cli_warn("The following regression jobs failed: {.val {failed}}.")
            for (i in failed) {
                cli_alert_warning("Job {.val {i}} failed. Rerunning locally")
                sge_res[[i]]$retv <- eval(parse(text = cmds[i]))
            }
        }
        cluster_models <- purrr::map(sge_res, "retv")
    } else {
        cluster_models <- safe_llply(seq_len(ncol(cluster_mat)), function(i) {
            cli_h1("Cluster {.val {i}}")
            regression_func(sequences, cluster_mat[, i], match_with_db = match_with_db, parallel = FALSE, ...)
        }, .parallel = parallel)
    }

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

    multiple_motifs <- "models" %in% names(cluster_models[[1]])

    if (!multiple_motifs) {
        motif_dataset <- purrr::imap_dfr(cluster_models, ~ .x$pssm %>% mutate(motif = .y)) %>%
            select(motif, pos, A, C, G, T)
        spat_dataset <- purrr::imap_dfr(cluster_models, ~ .x$spat %>% mutate(motif = .y)) %>%
            select(motif, bin, spat_factor)
    } else {
        motif_dataset <- purrr::imap_dfr(cluster_models, function(m, motif) {
            purrr::imap_dfr(m$models, ~ .x$pssm %>% mutate(cluster = motif, model = .y, motif = paste0(motif, ".", .y)))
        })
        spat_dataset <- purrr::imap_dfr(cluster_models, function(m, motif) {
            purrr::imap_dfr(m$models, ~ .x$spat %>% mutate(cluster = motif, model = .y, motif = paste0(motif, ".", .y)))
        })
    }

    res <- list(
        models = cluster_models,
        cluster_mat = cluster_mat,
        pred_mat = pred_mat,
        motif_dataset = motif_dataset,
        spat_dataset = spat_dataset,
        stats = stats
    )

    if (match_with_db) {
        if (!multiple_motifs) {
            cli_alert_info("Matching with motif databases")
            res$stats$db_match <- purrr::map_chr(cluster_models, "db_match")
            res$stats$db_match_cor <- purrr::map_dbl(cluster_models, "db_match_cor")
            res$pred_mat_db <- purrr::map(cluster_models, "db_match_pred") %>% do.call(cbind, .)
            colnames(res$pred_mat_db) <- names(res$stats$db_match)
            rownames(res$pred_mat_db) <- names(sequences)
            res$db_dataset <- purrr::imap_dfr(cluster_models, ~ .x$db_match_pssm %>% mutate(cluster = .y, motif = .x$db_match)) %>%
                select(cluster, motif, pos, A, C, G, T)
        } else {
            cli_alert_info("Cannot match with motif databases when {.field motif_num} > 1")
        }
    }

    if (multiple_motifs) {
        res$multi_stats <- purrr::imap_dfr(cluster_models, ~ .x$multi_stats %>% mutate(cluster = .y)) %>%
            select(cluster, everything())
    }

    if (screen_db) {
        cli_alert_info("Screening motif databases for {.val {ncol(cluster_mat)}} clusters")
        db_match <- screen_pwm.clusters(sequences, clusters, min_D = min_D, dataset = dataset, motifs = motifs, parallel = parallel, prior = prior, alternative = alternative)
        db_match <- purrr::imap_dfr(db_match, ~ tibble(cluster = .y, ks_D_db = max(.x), db_motif = rownames(db_match)[which.max(.x)]))
        res$stats <- res$stats %>% left_join(db_match, by = "cluster")
        for (clust in names(cluster_models)) {
            motif <- res$stats$db_motif[res$stats$cluster == clust]
            res$models[[clust]]$db_motif <- motif
            res$models[[clust]]$db_motif_score <- res$stats$ks_D_db[res$stats$cluster == clust]
            res$models[[clust]]$db_motif_pssm <- get_motif_pssm(motif)
            res$models[[clust]]$db_motif_pred <- extract_pwm(sequences, motif, prior = prior)[, 1]
        }
    }

    return(res)
}

#' Screen for motifs in a database for every cluster
#'
#' @param sequences a vector with the sequences
#' @param clusters a vector with the cluster assignments
#' @param min_D minimum distance to consider a match
#' @param only_best if TRUE, only return the best match for each cluster
#' @param alternative alternative hypothesis for the KS test. Can be "two.sided", "less" or "greater"
#'
#' @return a matrix with the KS D statistics for each cluster (columns) and every motif (rows)
#' that had at least one cluster with D >= min_D. If \code{only_best} is TRUE, a named vector
#' with the name of best motif match for each cluster is returned (regardless of \code{min_D}).
#'
#' @examples
#' \dontrun{
#' D_mat <- screen_pwm.clusters(cluster_sequences_example, clusters_example)
#' dim(D_mat)
#' D_mat[1:5, 1:5]
#'
#' # return only the best match
#' screen_pwm.clusters(cluster_sequences_example, clusters_example, only_best = TRUE)
#' }
#'
#' @inheritParams extract_pwm
#' @inheritDotParams compute_pwm
#' @export
screen_pwm.clusters <- function(sequences, clusters, dataset = all_motif_datasets(), motifs = NULL, parallel = getOption("prego.parallel", TRUE), min_D = 0.4, only_best = FALSE, prior = 0.01, alternative = "two.sided", ...) {
    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    sequences <- toupper(sequences)
    cluster_ids <- unique(clusters)

    res <- safe_daply(dataset, "motif", function(x) {
        pwm <- compute_pwm(sequences, x, prior = prior, ...)
        purrr::map_dbl(cluster_ids, function(cl) {
            suppressWarnings(ks.test(pwm[clusters == cl], pwm[clusters != cl], alternative = alternative)$statistic)
        })
    }, .parallel = parallel)
    if (is.vector(res)) {
        res <- matrix(res, ncol = length(res), nrow = 1)
    }
    colnames(res) <- cluster_ids
    res <- as.matrix(res)

    if (only_best) {
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

    return(as.data.frame(res))
}
