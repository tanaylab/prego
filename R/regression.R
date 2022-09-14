#' Perform a PWM regression
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper})
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param motif Initial motif to start the regression from. Can be either a string with a kmer where the character "*" indicates a
#' wildcard or a data frame with a pre-computed PSSM (see thre slot \code{pssm} in the return value of this function).
#' If NULL - a K-mer screen would be performed in order to find the best kmer for initialization.
#' @param motif_length Length of the seed motif. If the motif is shorter than this, it will be extended by wildcards (stars). Note that If the motif is longer than this, it will \emph{not} be truncated.
#' @param score_metric metric to use for optimizing the PWM. One of "r2" or "ks". When using "ks" the response variable should be a single vector of 0 and 1.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#' @param spat_min start of the spatial model from the beginning of the sequence (in bp)
#' @param spat_max end of the spatial model from the beginning of the sequence (in bp). If NULL - the spatial model
#' would end at the end of the sequence.
#' @param spat_bin size of the spatial bin (in bp).
#' @param improve_epsilon minimum improve in the objective function to continue the optimization
#' @param min_nuc_prob minimum nucleotide probability in every iteration
#' @param unif_prior uniform prior for nucleotide probabilities
#' @param include_response include the response in the resulting list (default: TRUE)
#' @param verbose show verbose messages.
#' @param seed random seed
#' @param consensus_single_thresh,consensus_double_thresh thresholds for the consensus sequence calculation
#' (single and double nucleotides)
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
#' }
#'
#' @examples
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
#' # binary response
#' res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' plot_regression_qc(res_binary)
#'
#' @inheritParams screen_kmers
#' @inheritDotParams screen_kmers
#' @export
regress_pwm <- function(sequences,
                        response,
                        motif = NULL,
                        motif_length = 15,
                        score_metric = "r2",
                        bidirect = TRUE,
                        spat_min = 0,
                        spat_max = NULL,
                        spat_bin = 50,
                        improve_epsilon = 0.0001,
                        min_nuc_prob = 0.001,
                        unif_prior = 0.05,
                        is_train = NULL,
                        include_response = TRUE,
                        seed = 60427,
                        verbose = FALSE,
                        kmer_length = 8,
                        motif_num = 1,
                        consensus_single_thresh = 0.6,
                        consensus_double_thresh = 0.85,
                        ...) {
    if (motif_num > 1) {
        return(regress_multiple_motifs(
            sequences = sequences,
            response = response,
            motif = motif,
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
            include_response = include_response,
            seed = seed,
            verbose = verbose,
            kmer_length = kmer_length,
            motif_num = motif_num,
            ...
        ))
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

    n_in_train <- sum(is_train)

    if (is.null(spat_max)) {
        spat_max <- nchar(sequences[1])
    }

    if (score_metric == "ks") {
        if (!is_binary_response(response)) {
            cli_abort("When {.field score_metric} is {.val ks}, {.field response} should be a single vector of 0 and 1")
        }
    }

    cli_alert_info("Number of response variables: {.val {ncol(response)}}")

    # get motif for initialization (either kmers or pssm)
    kmers <- NULL
    if (is.data.frame(motif)) { # initiazlie with pre-computed PSSM
        pssm <- motif
        if (!all(c("pos", "A", "C", "G", "T") %in% colnames(pssm))) {
            cli_abort("The {.field motif} PSSM data frame should have columns {.val pos}, {.val A}, {.val C}, {.val G}, {.val T}")
        }

        pssm <- pssm %>%
            arrange(as.numeric(pos)) %>%
            as.data.frame() %>%
            tibble::column_to_rownames("pos") %>%
            select(A, C, G, T) %>%
            as.matrix()

        motif <- ""
        cli_alert_info("Initializing regression with pre-computed PSSM")
    } else { # initialize with kmer
        pssm <- matrix()
        if (is.null(motif)) {
            cli_alert_info("Screening for kmers in order to initialize regression")
            kmers <- screen_kmers(sequences, response, kmer_length = kmer_length, ...)
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
        res$ks <- suppressWarnings(ks.test(res$pred[as.logical(response[, 1])], res$pred[!as.logical(response[, 1])]))
    }

    if (!is.null(kmers)) {
        res$kmers <- kmers
    }

    if (is_binary_response(response)) {
        cli_alert_success("Finished running regression. KS test D: {.val {round(res$ks$statistic, digits=2)}}, p-value: {.val {res$ks$p.value}}")
    } else {
        cli_alert_success("Finished running regression. R^2: {.val {round(res$r2, digits=4)}}")
    }

    return(res)
}
