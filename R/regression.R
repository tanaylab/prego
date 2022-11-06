#' Perform a PWM regression
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper})
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param motif Initial motif to start the regression from. Can be either a string with a kmer where the character "*" indicates a
#' wildcard or a data frame with a pre-computed PSSM (see the slot \code{pssm} in the return value of this function).
#' If NULL - a K-mer screen would be performed in order to find the best kmer for initialization. If \code{init_from_dataset} is TRUE, the regression would be initialized from the PSSM of the best motif in the dataset.
#' @param init_from_dataset initialize the regression from the PSSM of the best motif in \code{motif_dataset}, using \code{final_metric} as the metric. If TRUE, the \code{motif} parameter would be ignored. See \code{\link{screen_pwm}} for more details.
#' @param motif_length Length of the seed motif. If the motif is shorter than this, it will be extended by wildcards (stars). Note that If the motif is longer than this, it will \emph{not} be truncated.
#' @param score_metric metric to use for optimizing the PWM. One of "r2" or "ks". When using "ks" the response variable should be a single vector of 0 and 1.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#' @param spat_min start of the spatial model from the beginning of the sequence (in bp)
#' @param spat_max end of the spatial model from the beginning of the sequence (in bp). If NULL - the spatial model
#' would end at the end of the sequence.
#' @param spat_bin size of the spatial bin (in bp).
#' @param spat_model a previously computed spatial model (see \code{spat}) in the return value of this function. This can only be used when \code{motif} is a previously computed PSSM.
#' @param improve_epsilon minimum improve in the objective function to continue the optimization
#' @param min_nuc_prob minimum nucleotide probability in every iteration
#' @param unif_prior uniform prior for nucleotide probabilities
#' @param include_response include the response in the resulting list (default: TRUE)
#' @param verbose show verbose messages.
#' @param seed random seed
#' @param consensus_single_thresh,consensus_double_thresh thresholds for the consensus sequence calculation
#' (single and double nucleotides)
#' @param match_with_db match the resulting PWMs with motif databases using \code{pssm_match}. Note that the closest match
#' is returned, even if it is not similar enough in absolute terms.
#' @param motif_dataset  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs}, \code{JASPAR_motifs} or all_motif_datasets(). By default all_motif_datasets() would be used.
#' @param multi_kmers if TRUE, different candidates of kmers would be regressed in order to find the best seed according to \code{final_metric}.
#' @param final_metric metric to use in order to choose the best motif. One of 'ks' or 'r2'. Note that unlike \code{score_metric} which is used in the regression itself, this metric is used only for choosing the best motif out of all the runs on the sampled dataset.
#' @param kmer_length a vector of kmer lengths to screen in order to find the best seed motif.
#' @param max_cands maximum number of kmer candidates to try.
#' @param parallel whether to run optimization in parallel. use \code{set_parallel}
#' to set the number of cores to use.
#' @param motif_num Number of motifs to infer. When \code{motif_num} > 1, the function would run \code{motif_num} times, each time on the residuals of a linear model of all the previous runs (see \code{smooth_k} parameter). The best motif is then returned, while all the others are stored at 'models' in the return value.
#' @param smooth_k k for smoothing the predictions of each model in order to compute the residuals when \code{motif_num} > 1. The residuals are computed as \code{response} - running mean of size 'k' of the current model.
#' @param min_kmer_cor minimal correlation between the kmer and the response in order to use it as a seed.
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{pssm: }{data frame with the pssm matrix with the inferred motif, where rows are positions and columns are nucleotides.}
#' \item{spat: }{a data frame with the inferred spatial model, with the spatial factor for each bin. The bins are defined such that the first bin starts at \code{spat_min} and the last bin ends at \code{spat_max}, with a bin size of \code{spat_bin}.}
#' \item{pred: }{a vector with the predicted pwm for each sequence.}
#' \item{consensus: }{Consensus sequence based on the PSSM.}
#' \item{response: }{The response matrix. If \code{include_response} is FALSE, the response matrix is not included in the list.}
#' \item{r2: }{\eqn{r^2} of the prediction with respect to the each response variable.}
#' \item{ks: }{If response is binary, Kolmogorov-Smirnov test results of the predictions where the response was 1 vs the predictions where the response was 0.}
#' \item{seed_motif: }{The seed motif that started the regression.}
#' \item{kmers: }{The k-mers that were screened in order to find the best seed motif (if motif was NULL).}
#' \item{sample_idxs: }{The indices of the sequences that were used for the regression (only for \code{regress_pwm.sample}).}
#' }
#'
#' When \code{match_with_db} is TRUE, the following additional elements are returned:
#' \itemize{
#' \item{motif_db: }{The motif database that the most similar to the resulting PSSM.}
#' \item{db_match_dist: }{The distance between the resulting PSSM and the closest match in the motif database (KL divergence).}
#' \item{db_match_pssm: }{The PSSM of the closest match in the motif database.}
#' \item{db_match_pred: }{The predicted PWM of the closest match in the motif database.}
#' \item{db_match_r2: }{The \eqn{r^2} of the predicted PWM of the closest match in the motif database and the response}
#' \item{db_match_ks: }{If response is binary, the Kolmogorov-Smirnov test results of the predicted PWM of the closest match in the motif database where the response was 1 vs the predictions where the response was 0.}
#' }
#'
#' When \code{n_motifs} is greater than 1, the following additional elements are returned:
#' \itemize{
#' \item{models: }{A list (as above) of each inferred model}
#' \item{multi_stats: }{A data frame with the following columns: \code{model}, \code{score} (KS for binary, r^2 otherwise), \code{comb_score} (score for the combined linear model for models 1:i) and additional statistics per model}
#' }
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(sequences_example, response_mat_example)
#' res$pssm
#' res$spat
#' head(res$pred)
#'
#' plot_regression_qc(res)
#'
#' # intialize with a pre-computed PSSM
#' res1 <- regress_pwm(sequences_example, response_mat_example, motif = res$pssm)
#'
#' # intialize with a pre-computed PSSM and spatial model
#' res2 <- regress_pwm(
#'     sequences_example,
#'     response_mat_example,
#'     motif = res$pssm,
#'     spat_model = res$spat
#' )
#'
#' # binary response
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' plot_regression_qc(res_binary)
#'
#' # match with db
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], match_with_db = TRUE)
#' plot_regression_qc(res_binary)
#'
#' # initialize with a motif from the database
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], init_from_dataset = TRUE)
#'
#' # use multiple kmer seeds
#' res_multi <- regress_pwm(
#'     cluster_sequences_example,
#'     cluster_mat_example[, 1],
#'     multi_kmers = TRUE,
#'     kmer_length = 6:8,
#'     final_metric = "ks"
#' )
#' plot_regression_qc(res_multi)
#'
#' # Screen for multiple motifs
#' res_multi <- regress_pwm(
#'     cluster_sequences_example,
#'     cluster_mat_example[, 1],
#'     motif_num = 3,
#'     match_with_db = TRUE
#' )
#' res_multi$multi_stats
#' plot_regression_qc_multi(res_multi)
#' }
#'
#' # regress_multiple_motifs is an alias for regress_pwm with motif_num > 1
#' res_multi2 <- regress_multiple_motifs(
#'     cluster_sequences_example,
#'     cluster_mat_example[, 1],
#'     motif_num = 5,
#'     match_with_db = TRUE
#' )
#'
#' @inheritParams screen_kmers
#' @inheritDotParams screen_kmers
#' @export
regress_pwm <- function(sequences,
                        response,
                        motif = NULL,
                        motif_length = 15,
                        init_from_dataset = FALSE,
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
                        final_metric = "r2",
                        max_cands = 10,
                        min_gap = 0,
                        max_gap = 1,
                        min_kmer_cor = 0.1,
                        motif_num = 1,
                        smooth_k = 100,
                        consensus_single_thresh = 0.6,
                        consensus_double_thresh = 0.85,
                        match_with_db = FALSE,
                        motif_dataset = all_motif_datasets(),
                        parallel = getOption("prego.parallel", FALSE),
                        ...) {
    if (motif_num > 1) {
        return(
            regress_multiple_motifs(
                sequences = sequences,
                response = response,
                motif = motif,
                motif_length = motif_length,
                score_metric = score_metric,
                bidirect = bidirect,
                spat_min = spat_min,
                spat_max = spat_max,
                spat_bin = spat_bin,
                spat_model = spat_model,
                improve_epsilon = improve_epsilon,
                min_nuc_prob = min_nuc_prob,
                unif_prior = unif_prior,
                is_train = is_train,
                include_response = include_response,
                seed = seed,
                verbose = verbose,
                kmer_length = kmer_length,
                multi_kmers = multi_kmers,
                final_metric = final_metric,
                max_cands = max_cands,
                min_gap = min_gap,
                max_gap = max_gap,
                min_kmer_cor = min_kmer_cor,
                motif_num = motif_num,
                smooth_k = smooth_k,
                consensus_single_thresh = consensus_single_thresh,
                consensus_double_thresh = consensus_double_thresh,
                match_with_db = match_with_db,
                motif_dataset = motif_dataset,
                parallel = parallel,
                ...
            )
        )
    }

    if (!(score_metric %in% c("r2", "ks"))) {
        cli_abort("score_metric must be one of {.val r2} or {.val ks}")
    }

    if (is.null(is_train)) {
        is_train <- rep(TRUE, length(sequences))
    }

    if (length(sequences) != length(is_train)) {
        cli_abort("The number of sequences and the length of {.field is_train} vector do not match")
    }

    if (is.null(nrow(response))) {
        response <- matrix(response, ncol = 1)
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    n_in_train <- sum(is_train)

    if (is.null(spat_max)) {
        spat_max <- nchar(sequences[1])
    }

    if (spat_min < 0) {
        cli_abort("{.field spat_min} must be non-negative")
    }

    if ((spat_max > nchar(sequences[1])) | (spat_max < spat_min)) {
        cli_abort("{.field spat_max} must be between {.field spat_min} and the length of the sequences")
    }

    if (as.integer((spat_max - spat_min) / spat_bin) != (spat_max - spat_min) / spat_bin) {
        cli_abort("{.field spat_bin} must be a divisor of {.field spat_max} - {.field spat_min}")
    }

    if (!is.null(spat_model)) {
        if (!is.data.frame(motif)) {
            cli_abort("If {.field spat_model} is provided, {.field motif} must be a previously computed PSSM")
        }
        validate_spat(spat_model)
        spat_bin <- unique(diff(spat_model$bin))
        spat_model <- spat_model$spat_factor
    }

    if (score_metric == "ks") {
        if (!is_binary_response(response)) {
            cli_abort("When {.field score_metric} is {.val ks}, {.field response} should be a single vector of 0 and 1")
        }
    }

    cli_alert_info("Number of response variables: {.val {ncol(response)}}")

    if (init_from_dataset) {
        cli_alert_info("Initializing from dataset")
        motif_name <- screen_pwm(sequences, response, metric = final_metric, prior = unif_prior, bidirect = bidirect, only_best = TRUE)
        motif <- get_motif_pssm(motif_name$motif)
        motif <- pssm_add_prior(motif, prior = unif_prior)
        cli_alert_info("Best motif from dataset: {.val {motif_name$motif}}")
    }

    if (!is.null(motif) && multi_kmers) {
        cli_warn("Motif is provided, {.field multi_kmers} will be ignored")
    }

    # get motif for initialization (either kmers or pssm)
    kmers <- NULL
    if (is.data.frame(motif)) { # initiazlie with pre-computed PSSM
        pssm <- motif
        if (!all(c("pos", "A", "C", "G", "T") %in% colnames(pssm))) {
            cli_abort("The {.field motif} PSSM data frame should have columns {.val pos}, {.val A}, {.val C}, {.val G}, {.val T}")
        }

        pssm <- pssm_to_mat(pssm)
        consensus <- consensus_from_pssm(pssm, consensus_single_thresh, consensus_double_thresh)

        motif <- ""
        if (!is.na(consensus)) {
            cli_alert_info("Initializing regression with pre-computed PSSM, consensus: {.val {consensus_from_pssm(pssm)}}")
        } else {
            cli_alert_info("Initializing regression with pre-computed PSSM")
        }
    } else { # initialize with kmer
        pssm <- matrix()
        if (is.null(motif)) {
            if (multi_kmers) {
                return(regress_pwm.multi_kmers(
                    sequences = sequences,
                    response = response,
                    motif_length = motif_length,
                    score_metric = score_metric,
                    bidirect = bidirect,
                    spat_min = spat_min,
                    spat_max = spat_max,
                    spat_bin = spat_bin,
                    spat_model = spat_model,
                    improve_epsilon = improve_epsilon,
                    min_nuc_prob = min_nuc_prob,
                    unif_prior = unif_prior,
                    is_train = is_train,
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
                    final_metric = final_metric,
                    parallel = parallel,
                    match_with_db = match_with_db,
                    motif_dataset = motif_dataset,
                    ...
                ))
            }
            cli_alert_info("Screening for kmers in order to initialize regression")
            kmers <- screen_kmers(sequences, response, kmer_length = kmer_length, min_gap = min_gap, max_gap = max_gap, min_cor = min_kmer_cor, ...)
            motif <- kmers$kmer[which.max(abs(kmers$max_r2))]
            if (length(motif) == 0) { # could not find any kmer
                motif <- paste(rep("*", motif_length), collapse = "")
                cli_alert_info("Could not find any kmer. Initializing with {.val {motif}}")
            }
        }
        if (stringr::str_length(motif) < motif_length) {
            motif <- stringr::str_pad(motif, motif_length, "*", side = "both")
            cli_alert_info("Motif is shorter than {.val {motif_length}}, extending with wildcards")
        }

        cli_alert_info("Initializing regression with the following motif: {.val {motif}}")
    }

    cli_alert_info("Running regression")
    cli_ul(c(
        "Motif length: {.val {motif_length}}",
        "Bidirectional: {.val {bidirect}}",
        "Spat min: {.val {spat_min}}",
        "Spat max: {.val {spat_max}}",
        "Spat bin: {.val {spat_bin}}",
        "Improve epsilon: {.val {improve_epsilon}}",
        "Min nuc prob: {.val {min_nuc_prob}}",
        "Uniform prior: {.val {unif_prior}}",
        "Score metric: {.val {score_metric}}",
        "Seed: {.val {seed}}"
    ))

    res <- regress_pwm_cpp(
        toupper(sequences),
        response,
        is_train,
        motif = motif,
        spat_min = spat_min,
        spat_max = spat_max,
        min_nuc_prob = min_nuc_prob,
        spat_bin = spat_bin,
        spat_factor = spat_model,
        improve_epsilon = improve_epsilon,
        is_bidirect = bidirect,
        unif_prior = unif_prior,
        score_metric = score_metric,
        verbose = verbose,
        seed = seed,
        pssm_mat = pssm,
        consensus_single_thresh = consensus_single_thresh,
        consensus_double_thresh = consensus_double_thresh
    )


    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2
    res$seed_motif <- motif

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])], alternative = "less"))
    }

    if (!is.null(kmers)) {
        res$kmers <- kmers
    }

    res$consensus <- consensus_from_pssm(res$pssm, consensus_single_thresh, consensus_double_thresh)
    cli_alert_success("Finished running regression. Consensus: {.val {res$consensus}}")



    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset)
    }

    if (is_binary_response(response)) {
        cli_alert_success("KS test D: {.val {round(res$ks$statistic, digits=4)}}, p-value: {.val {res$ks$p.value}}")
    } else {
        cli_alert_success("R^2: {.val {round(res$r2, digits=4)}}")
    }

    return(res)
}

add_regression_db_match <- function(reg, sequences, motif_dataset, parallel = getOption("prego.parallel", TRUE)) {
    best_match <- pssm_match(reg$pssm, motif_dataset, parallel = parallel)[1, ]
    reg$db_match <- best_match$motif
    reg$db_match_dist <- best_match$dist
    cli_alert_info("Best match in the database: {.val {best_match$motif}}, KL: {.val {round(best_match$dist, digits = 3)}}")
    reg$db_match_pssm <- motif_dataset %>%
        filter(motif == best_match$motif) %>%
        select(pos, A, C, G, T)
    reg$db_match_pred <- compute_pwm(sequences, reg$db_match_pssm)
    reg$db_match_r2 <- tgs_cor(reg$response, as.matrix(reg$db_match_pred))[, 1]^2
    if (is_binary_response(reg$response)) {
        reg$db_match_ks <- suppressWarnings(ks.test(reg$db_match_pred[as.logical(reg$response[, 1])], reg$db_match_pred[!as.logical(reg$response[, 1])], alternative = "less"))
        cli_alert_success("{.val {reg$db_match}} KS test D: {.val {round(reg$db_match_ks$statistic, digits=4)}}, p-value: {.val {reg$db_match_ks$p.value}}")
    }
    return(reg)
}

regress_pwm.multi_kmers <- function(sequences,
                                    response,
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
                                    kmer_length = 6:8,
                                    max_cands = 10,
                                    min_gap = 0,
                                    max_gap = 1,
                                    min_kmer_cor = 0.1,
                                    consensus_single_thresh = 0.6,
                                    consensus_double_thresh = 0.85,
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

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    regress_pwm_single_kmer <- purrr::partial(
        regress_pwm,
        sequences = sequences,
        response = response,
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
        consensus_single_thresh = consensus_single_thresh,
        consensus_double_thresh = consensus_double_thresh,
        match_with_db = FALSE
    )

    cli_h3("Generate candidate kmers")
    cand_kmers <- get_cand_kmers(sequences, response, kmer_length, min_gap, max_gap, min_kmer_cor, verbose, parallel, max_cands = max_cands, ...)

    cli_h3("Regress each candidate kmer on sampled data")
    cli_alert_info("Running regression on {.val {length(cand_kmers)}} candidate kmers")
    cli_ul(c(
        "Bidirectional: {.val {bidirect}}",
        "Spat min: {.val {spat_min}}",
        "Spat max: {.val {spat_max}}",
        "Spat bin: {.val {spat_bin}}",
        "Improve epsilon: {.val {improve_epsilon}}",
        "Min nuc prob: {.val {min_nuc_prob}}",
        "Uniform prior: {.val {unif_prior}}",
        "Score metric: {.val {score_metric}}",
        "Seed: {.val {seed}}"
    ))
    res_kmer_list <- plyr::llply(cli_progress_along(cand_kmers), function(i) {
        motif <- cand_kmers[i]
        cli_alert("regressing with seed: {.val {motif}}")
        r <- regress_pwm_single_kmer(motif = motif) %>%
            suppressMessages()
        if (final_metric == "ks") {
            if (!is_binary_response(response)) {
                cli_abort("Cannot use {.field final_metric} {.val ks} when {.field response} is not binary")
            }
            r$score <- r[[final_metric]]$statistic
        } else if (final_metric == "r2") {
            r$score <- r[[final_metric]]
        } else {
            cli_abort("Unknown {.field final_metric} (can be 'ks' or 'r2')")
        }
        cli_alert("{.val {motif}}, score ({final_metric}): {.val {r$score}}")
        return(r)
    }, .parallel = parallel)

    scores <- sapply(res_kmer_list, function(x) x$score)

    if (length(which.max(scores)) == 0) {
        cli_alert_warning("No motifs found")
        res <- res_kmer_list[[1]]
    } else {
        res <- res_kmer_list[[which.max(scores)]]
    }

    res$kmers <- cand_kmers

    if (include_response) {
        res$response <- response
    }

    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset, parallel = parallel)
    }

    cli_alert_info("Best motif: {.val {res$seed_motif}}, score ({final_metric}): {.val {max(scores)}}")

    return(res)
}

get_cand_kmers <- function(sequences, response, kmer_length, min_gap, max_gap, min_kmer_cor, verbose, parallel = FALSE, max_cands = 10, ...) {
    all_kmers <- plyr::ldply(cli_progress_along(kmer_length), function(i) {
        screen_kmers(sequences, response, kmer_length = kmer_length[i], min_gap = min_gap, max_gap = max_gap, ...) %>%
            mutate(len = kmer_length[i], verbose = FALSE) %>%
            suppressMessages()
    }, .parallel = parallel)

    best_kmer <- all_kmers$kmer[which.max(abs(all_kmers$max_r2))] # return at least one kmer

    if (length(best_kmer) == 0) { # could not find any kmer
        res <- paste(rep("*", kmer_length), collapse = "")
        cli_alert_info("Could not find any kmer. Initializing with {.val {res}}")
        return(res)
    }

    all_kmers <- all_kmers %>%
        # filter by correlation
        filter(sqrt(max_r2) > min_kmer_cor) %>%
        dplyr::distinct(kmer, .keep_all = TRUE)

    cands <- all_kmers %>%
        slice_max(n = min(nrow(all_kmers), max_cands), order_by = abs(max_r2)) %>%
        arrange(desc(abs(max_r2)))


    dist_mat <- stringdist::stringdistmatrix(cands$kmer, cands$kmer, method = "osa", nthread = 1)
    dist_mat[dist_mat != 1] <- NA
    if (ncol(dist_mat) == nrow(dist_mat)) {
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
