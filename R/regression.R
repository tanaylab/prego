#' Perform a PWM regression
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper}). Please make sure that the sequences are long enough to cover \code{spat_num_bins} * \code{spat_bin_size} bp, and that they are centered around the motif/signal.
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param motif Initial motif to start the regression from. Can be either a string with a kmer where the character "*" indicates a
#' wildcard or a data frame with a pre-computed PSSM (see the slot \code{pssm} in the return value of this function).
#' If NULL - a K-mer screen would be performed in order to find the best kmer for initialization. If \code{init_from_dataset} is TRUE, the regression would be initialized from the PSSM of the best motif in the dataset.
#' @param init_from_dataset initialize the regression from the PSSM of the best motif in \code{motif_dataset}, using \code{final_metric} as the metric. If TRUE, the \code{motif} parameter would be ignored. See \code{\link{screen_pwm}} for more details.
#' @param motif_length Length of the seed motif. If the motif is shorter than this, it will be extended by wildcards (stars). Note that If the motif is longer than this, it will \emph{not} be truncated.
#' @param score_metric metric to use for optimizing the PWM. One of "r2" or "ks". When using "ks" the response variable should be a single vector of 0 and 1.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#' @param spat_bin_size size of the spatial bin (in bp).
#' @param spat_num_bins number of spatial bins. Please make sure that the sequences are long enough to cover this number of bins. bp outside of spat_bin_size * spat_num_bins would be ignored. If \code{bidirect} is TRUE, the number of bins should be odd as 'prego' symmetrizes the motif around the center bin.
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
#' @param screen_db Screen \code{motif_dataset} using \code{screen_pwm} and use the best motif as the initial motif. If TRUE, the following fields would be added to the return value:
#' "db_motif", "db_motif_pred", "db_motif_pssm" and "db_motif_score".
#' @param motif_dataset  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs}, \code{JASPAR_motifs} or all_motif_datasets(). By default all_motif_datasets() would be used.
#' @param multi_kmers if TRUE, different candidates of kmers would be regressed in order to find the best seed according to \code{final_metric}.
#' @param final_metric metric to use in order to choose the best motif. One of 'ks' or 'r2'. Note that unlike \code{score_metric} which is used in the regression itself, this metric is used only for choosing the best motif out of all the runs on the sampled dataset. If NULL - 'ks' would be used for binary response and 'r2' for continuous response.
#' @param kmer_length a vector of kmer lengths to screen in order to find the best seed motif.
#' @param max_cands maximum number of kmer candidates to try.
#' @param parallel whether to run optimization in parallel. use \code{set_parallel}
#' to set the number of cores to use.
#' @param motif_num Number of motifs to infer. When \code{motif_num} > 1, the function would run \code{motif_num} times, each time on the residuals of a linear model of all the previous runs (see \code{smooth_k} parameter). The best motif is then returned, while all the others are stored at 'models' in the return value.
#' @param smooth_k k for smoothing the predictions of each model in order to compute the residuals when \code{motif_num} > 1. The residuals are computed as \code{response} - running mean of size 'k' of the current model.
#' @param min_kmer_cor minimal correlation between the kmer and the response in order to use it as a seed.
#' @param internal_num_folds number of folds to use in the internal cross-validation.
#' @param alternative alternative hypothesis for the p-value calculation when using \code{ks.test}. One of "two.sided", "less" or "greater".
#' @param sample_for_kmers Use a random sample of the dataset in order to find the best kmer. This is useful when the dataset is very large and the kmer screen would take a long time. Note that the final regression would be performed on the entire dataset. Only relevant when \code{multi_kmers} is TRUE.
#' @param sample_frac fraction of the dataset to use for the kmer screen. Default: 0.1.
#' @param sample_idxs indices of the sequences to use for the kmer screen. If NULL, a random sample would be used.
#' @param sample_ratio ratio between the '1' category and the '0' category in the sampled dataset (for binary response). Relevant only when \code{sample_frac} is NULL.
#' @param log_energy transform the energy to log scale on each iteration.
#' @param energy_func a function to transform the energy at each iteration. Should accept a numeric vector and return a numeric vector. e.g. \code{log} or \code{function(x) x^2}. Note that the range of the input energies is between 0 and 1 (the probability of the motif in the sequence), so if you inferred the function using the the returned energies (which are in log scale) you should make sure that the function first log transforms using \code{log_energy=TRUE}.
#' @param xmin,xmax,npts range for the energy function and the number of points to use for its interpolation.
#' @param energy_func_generator a function to generate the energy function when regressing multiple motifs. Should accept the result of the previous iteration + the original response and return a function similar to \code{energy_func}. e.g. \code{
#' function(prev_reg, resp) {
#'        df <- data.frame(x = prev_reg$pred, y = resp)
#'        fn_gam <- as.formula("y ~ s(x, k=3, bs='cr')")
#'        model <- mgcv::gam(fn_gam, family = binomial(link = "logit"), data = df, method="REML")
#'        function(z){
#'            mgcv::predict.gam(object = model, newdata = data.frame(x = z))
#' }}}.
#' When this parameter is not NULL, energy_func_generator would create an energy function and then run another step of regression initialized with the previous motif with \code{energy_func} as the energy function. This is useful when the energy function is not monotonic, for example - one might want to use a gam model to fit the energy function like in the example above.
#' @param optimize_pwm optimize the PWM model (Default: TRUE). If FALSE, the PWM model would be used as the initial model for the spatial model.
#' @param optimize_spat optimize the spatial model (Default: TRUE). If FALSE, the spatial model would be used as the initial model for the PWM model.
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{pssm: }{data frame with the pssm matrix with the inferred motif, where rows are positions and columns are nucleotides.}
#' \item{spat: }{a data frame with the inferred spatial model, with the spatial factor for each bin.}
#' \item{pred: }{a vector with the predicted pwm for each sequence.}
#' \item{consensus: }{Consensus sequence based on the PSSM.}
#' \item{response: }{The response matrix. If \code{include_response} is FALSE, the response matrix is not included in the list.}
#' \item{r2: }{\eqn{r^2} of the prediction with respect to the each response variable.}
#' \item{ks: }{If response is binary, Kolmogorov-Smirnov test results of the predictions where the response was 1 vs the predictions where the response was 0.}
#' \item{seed_motif: }{The seed motif that started the regression.}
#' \item{kmers: }{The k-mers that were screened in order to find the best seed motif (if motif was NULL).}
#' \item{sample_idxs: }{The indices of the sequences that were used for the regression (only for \code{regress_pwm.sample}).}
#' \item{predict: }{a function that can be used to predict the PWM for a new sequence.}
#' }
#'
#' When \code{match_with_db} is TRUE, the following additional elements are returned:
#' \itemize{
#' \item{motif_db: }{The motif database that the most similar to the resulting PSSM.}
#' \item{db_match_cor: }{The correlation between the resulting PSSM and the closest match in the motif database.}
#' \item{db_match_pssm: }{The PSSM of the closest match in the motif database.}
#' \item{db_match_pred: }{The predicted PWM of the closest match in the motif database.}
#' \item{db_match_r2: }{The \eqn{r^2} of the predicted PWM of the closest match in the motif database and the response}
#' \item{db_match_ks: }{If response is binary, the Kolmogorov-Smirnov test results of the predicted PWM of the closest match in the motif database where the response was 1 vs the predictions where the response was 0.}
#' }
#'
#' When \code{screen_db} is TRUE, the following additional elements are returned:
#' \itemize{
#' \item{db_motif: }{The best motif from the motif database.}
#' \item{db_motif_pred: }{The predicted PWM of the best motif from the motif database.}
#' \item{db_motif_pssm: }{The PSSM of the best motif from the motif database.}
#' \item{db_motif_score: }{The score of the best motif from the motif database.}
#' }
#'
#' When \code{n_motifs} is greater than 1, a list with the following elements is returned:
#' \itemize{
#' \item{models: }{A list (as above) of each inferred model}
#' \item{multi_stats: }{A data frame with the following columns: \code{model}, \code{score} (KS for binary, r^2 otherwise), \code{comb_score} (score for the combined linear model for models 1:i) and additional statistics per model}
#' \item{pred: }{a vector with the predicted pwm for using a linear model of the combined scores.}
#' \item{comb_modle: }{a linear model of the combined scores.}
#' \item{predict: }{a function that can be used to predict the PWM for a new sequence.}
#' \item{predict_multi: }{a function that can be used to predict the PWM for the different models for a new sequence}
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
#' # Screen for best db motif
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], screen_db = TRUE)
#' plot_regression_qc(res_binary)
#'
#' # initialize with a motif from the database
#' res_binary <- regress_pwm(
#'     cluster_sequences_example,
#'     cluster_mat_example[, 1],
#'     init_from_dataset = TRUE
#' )
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
                        spat_bin_size = NULL,
                        spat_num_bins = 7,
                        spat_model = NULL,
                        improve_epsilon = 0.0001,
                        min_nuc_prob = 0.001,
                        unif_prior = 0.05,
                        include_response = TRUE,
                        seed = 60427,
                        verbose = FALSE,
                        kmer_length = 8,
                        multi_kmers = TRUE,
                        final_metric = NULL,
                        max_cands = 10,
                        min_gap = 0,
                        max_gap = 1,
                        min_kmer_cor = 0.08,
                        motif_num = 1,
                        smooth_k = 100,
                        consensus_single_thresh = 0.6,
                        consensus_double_thresh = 0.85,
                        internal_num_folds = 1,
                        match_with_db = TRUE,
                        screen_db = FALSE,
                        motif_dataset = all_motif_datasets(),
                        parallel = getOption("prego.parallel", FALSE),
                        alternative = "less",
                        sample_for_kmers = FALSE,
                        sample_frac = NULL,
                        sample_idxs = NULL,
                        sample_ratio = 1,
                        log_energy = FALSE,
                        energy_func = NULL,
                        xmin = -100,
                        xmax = 100,
                        npts = 1e4,
                        energy_func_generator = NULL,
                        optimize_pwm = TRUE,
                        optimize_spat = TRUE,
                        ...) {
    set.seed(seed)
    if (motif_num > 1) {
        return(
            regress_multiple_motifs(
                sequences = sequences,
                response = response,
                motif = motif,
                motif_length = motif_length,
                score_metric = score_metric,
                bidirect = bidirect,
                spat_bin_size = spat_bin_size,
                spat_num_bins = spat_num_bins,
                spat_model = spat_model,
                improve_epsilon = improve_epsilon,
                min_nuc_prob = min_nuc_prob,
                unif_prior = unif_prior,
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
                internal_num_folds = internal_num_folds,
                match_with_db = match_with_db,
                motif_dataset = motif_dataset,
                parallel = parallel,
                alternative = alternative,
                sample_for_kmers = sample_for_kmers,
                sample_frac = sample_frac,
                sample_idxs = sample_idxs,
                sample_ratio = sample_ratio,
                log_energy = log_energy,
                energy_func = energy_func,
                energy_func_generator = energy_func_generator,
                xmin = xmin,
                xmax = xmax,
                npts = npts,
                optimize_pwm = optimize_pwm,
                optimize_spat = optimize_spat,
                ...
            )
        )
    }

    if (!(score_metric %in% c("r2", "ks"))) {
        cli_abort("score_metric must be one of {.val r2} or {.val ks}")
    }

    if (is.null(nrow(response))) {
        response <- matrix(response, ncol = 1)
    }

    if (!all(is.numeric(response))) {
        cli_abort("{.field response} must be numeric")
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    max_seq_len <- nchar(sequences[1])

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

    if (is.null(final_metric)) {
        if (is_binary_response(response)) {
            final_metric <- "ks"
        } else {
            final_metric <- "r2"
        }
    }
    cli_alert_info("Using {.val {final_metric}} as the final metric")

    cli_alert_info("Number of response variables: {.val {ncol(response)}}")

    if (!is.null(energy_func) && !log_energy) {
        cli_warn("Energy function was provided, but {.field log_energy} is {.val FALSE}. Note that the energies during the regression are not log-transformed.")
    }

    spat <- calc_spat_min_max(spat_num_bins, max_seq_len, spat_bin_size)

    if (init_from_dataset) {
        cli_alert_info("Initializing from dataset")
        motif_name <- screen_pwm(sequences, response, metric = final_metric, prior = unif_prior, bidirect = bidirect, only_best = TRUE, alternative = alternative)
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
                    spat_bin_size = spat_bin_size,
                    spat_num_bins = spat_num_bins,
                    spat_model = spat_model,
                    improve_epsilon = improve_epsilon,
                    min_nuc_prob = min_nuc_prob,
                    unif_prior = unif_prior,
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
                    final_metric = final_metric,
                    parallel = parallel,
                    match_with_db = match_with_db,
                    screen_db = screen_db,
                    motif_dataset = motif_dataset,
                    alternative = alternative,
                    sample_for_kmers = sample_for_kmers,
                    sample_frac = sample_frac,
                    sample_idxs = sample_idxs,
                    sample_ratio = sample_ratio,
                    log_energy = log_energy,
                    energy_func = energy_func,
                    xmin = xmin,
                    xmax = xmax,
                    npts = npts,
                    optimize_pwm = optimize_pwm,
                    optimize_spat = optimize_spat,
                    ...
                ))
            }
            cli_alert_info("Screening for kmers in order to initialize regression")
            kmers <- screen_kmers(sequences, response, kmer_length = kmer_length, min_gap = min_gap, max_gap = max_gap, min_cor = min_kmer_cor, ...)
            motif <- kmers$kmer[which.max(abs(kmers$max_r2))]
            if (length(motif) == 0) { # could not find any kmer
                cli_alert_info("Could not find any kmer with correlation above {.val {min_kmer_cor}}. Trying with a threshold of {.val {min_kmer_cor / 2}}")
                kmers <- screen_kmers(sequences, response, kmer_length = kmer_length, min_gap = min_gap, max_gap = max_gap, min_cor = min_kmer_cor / 2, ...)
                motif <- kmers$kmer[which.max(abs(kmers$max_r2))]
                if (length(motif) == 0) {
                    motif <- paste(rep("*", motif_length), collapse = "")
                    cli_alert_info("Could not find any kmer. Initializing with {.val {motif}}")
                }
            }
        }
        if (stringr::str_length(motif) < motif_length) {
            motif <- stringr::str_pad(motif, motif_length, "*", side = "both")
            cli_alert_info("Motif is shorter than {.val {motif_length}}, extending with wildcards")

            # substitute 'N' with wildcards
            motif <- stringr::str_replace_all(motif, "N", "*")
        }

        cli_alert_info("Initializing regression with the following motif: {.val {motif}}")
    }



    cli_alert_info("Running regression")
    cli_ul(c(
        "Motif length: {.val {motif_length}}",
        "Bidirectional: {.val {bidirect}}",
        "Spat min: {.val {spat$spat_min}}",
        "Spat max: {.val {spat$spat_max}}",
        "Spat bin size: {.val {spat_bin_size}}",
        "Number of bins: {.val {spat_num_bins}}",
        "Improve epsilon: {.val {improve_epsilon}}",
        "Min nuc prob: {.val {min_nuc_prob}}",
        "Uniform prior: {.val {unif_prior}}",
        "Score metric: {.val {score_metric}}",
        "Seed: {.val {seed}}"
    ))

    sequences <- stringr::str_sub(sequences, start = spat$spat_min, end = spat$spat_max - 1)

    is_train <- rep(TRUE, length(sequences))

    res <- regress_pwm_cpp(
        toupper(sequences),
        response,
        is_train,
        motif = motif,
        spat_bin = spat_bin_size,
        spat_min = 0,
        spat_max = nchar(sequences[1]),
        min_nuc_prob = min_nuc_prob,
        spat_factor = spat_model,
        improve_epsilon = improve_epsilon,
        is_bidirect = bidirect,
        unif_prior = unif_prior,
        score_metric = score_metric,
        verbose = verbose,
        seed = seed,
        pssm_mat = pssm,
        consensus_single_thresh = consensus_single_thresh,
        consensus_double_thresh = consensus_double_thresh,
        num_folds = internal_num_folds,
        log_energy = log_energy,
        energy_func = energy_func,
        xmin = xmin,
        xmax = xmax,
        npts = npts,
        optimize_pwm = optimize_pwm,
        optimize_spat = optimize_spat
    )


    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2
    res$seed_motif <- motif

    if (is_binary_response(response)) {
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])], alternative = alternative))
        res$score <- res$ks$statistic
    } else {
        res$score <- res$r2
    }

    if (!is.null(kmers)) {
        res$kmers <- kmers
    }

    res$consensus <- consensus_from_pssm(res$pssm, consensus_single_thresh, consensus_double_thresh)
    cli_alert_success("Finished running regression. Consensus: {.val {res$consensus}}")



    if (match_with_db) {
        res <- add_regression_db_match(res, sequences, motif_dataset, alternative = alternative)
    }

    if (screen_db) {
        res <- add_regression_db_screen(res, response, sequences, motif_dataset, final_metric, prior = unif_prior, bidirect = bidirect, alternative = alternative, parallel = parallel)
    }

    if (is_binary_response(response)) {
        cli_alert_success("KS test D: {.val {round(res$ks$statistic, digits=4)}}, p-value: {.val {res$ks$p.value}}")
    } else {
        cli_alert_success("R^2: {.val {round(res$r2, digits=4)}}")
    }

    res$spat_min <- spat$spat_min
    res$spat_max <- spat$spat_max
    res$spat_bin_size <- spat_bin_size
    res$bidirect <- bidirect
    res$seq_length <- nchar(sequences[1])

    res <- add_predict_function(res, spat, bidirect, energy_func)

    return(res)
}

add_predict_function <- function(res, spat, bidirect, energy_func) {
    func_env <- new.env()
    func_env$pssm <- res$pssm
    func_env$spat <- res$spat
    func_env$energy_func <- energy_func
    func_env$bidirect <- bidirect
    func_env$spat_min <- spat$spat_min
    func_env$spat_max <- spat$spat_max

    res$predict <- function(x) {
        x <- stringr::str_sub(x, start = spat_min, end = spat_max - 1)
        e <- compute_pwm(x, pssm, spat = spat, bidirect = bidirect, spat_min = 0, spat_max = nchar(x)[1], prior = 0)
        if (!is.null(energy_func)) {
            e <- energy_func(e)
        }
        return(e)
    }

    environment(res$predict) <- func_env

    return(res)
}

add_regression_db_screen <- function(res, response, sequences, motif_dataset, metric, prior, bidirect, alternative, parallel = getOption("prego.parallel", TRUE)) {
    scr <- screen_pwm(sequences, response, motif_dataset, metric = metric, prior = prior, bidirect = bidirect, parallel = parallel, alternative = alternative, only_best = TRUE)

    cli_alert_info("Best db motif: {.val {scr$motif[1]}}")
    cli_alert_info("Best db motif score ({.val {metric}}): {.val {scr$score[1]}}")
    res$db_motif <- scr$motif[1]
    res$db_motif_pred <- extract_pwm(sequences, scr$motif[1], prior = prior, bidirect = bidirect)[, 1]
    res$db_motif_pssm <- get_motif_pssm(scr$motif[1])
    res$db_motif_score <- scr$score[1]

    return(res)
}

add_regression_db_match <- function(reg, sequences, motif_dataset, alternative, parallel = getOption("prego.parallel", TRUE)) {
    best_match <- pssm_match(reg$pssm, motif_dataset, parallel = parallel)[1, ]
    reg$db_match <- best_match$motif
    reg$db_match_cor <- best_match$cor
    cli_alert_info("Best match in the database: {.val {best_match$motif}}, cor: {.val {round(best_match$cor, digits = 3)}}")
    reg$db_match_pssm <- motif_dataset %>%
        filter(motif == best_match$motif) %>%
        select(pos, A, C, G, T)
    reg$db_match_pred <- compute_pwm(sequences, reg$db_match_pssm)
    reg$db_match_r2 <- tgs_cor(as.matrix(reg$response), as.matrix(reg$db_match_pred))[, 1]^2
    if (is_binary_response(reg$response)) {
        reg$db_match_ks <- suppressWarnings(ks.test(reg$db_match_pred[as.logical(reg$response[, 1])], reg$db_match_pred[!as.logical(reg$response[, 1])], alternative = alternative))
        cli_alert_success("{.val {reg$db_match}} KS test D: {.val {round(reg$db_match_ks$statistic, digits=4)}}, p-value: {.val {reg$db_match_ks$p.value}}")
    }
    return(reg)
}



calc_spat_min_max <- function(spat_num_bins, max_seq_len, spat_bin_size = NULL) {
    if (is.null(spat_bin_size)) {
        # make sure it is an even number
        spat_bin_size <- round(max_seq_len / spat_num_bins / 2) * 2
    }
    if (spat_bin_size %% 2 != 0) {
        cli_abort("The {.field spat_bin_size} must be an even number")
    }
    if (spat_num_bins %% 2 != 1) {
        cli_abort("The {.field spat_num_bins} must be an odd number")
    }
    if (spat_bin_size * spat_num_bins > max_seq_len) {
        cli_abort("The {.field spat_bin_size} ({.val {spat_bin_size}}) times the {.field spat_num_bins} ({.val {spat_num_bins}}) must be smaller than the maximum sequence length ({.val {max_seq_len}})")
    }

    center <- round(max_seq_len / 2)

    if (spat_num_bins == 1) {
        spat_min <- center - spat_bin_size / 2
        spat_max <- center + spat_bin_size / 2
    } else {
        # position one bin at the center, and then add bins to the left and to the right
        spat_min <- center - ((spat_num_bins - 1) / 2) * spat_bin_size - spat_bin_size / 2
        spat_max <- center + ((spat_num_bins - 1) / 2) * spat_bin_size + spat_bin_size / 2
    }


    return(list(spat_min = round(spat_min), spat_max = round(spat_max)))
}
