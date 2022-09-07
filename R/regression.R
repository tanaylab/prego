#' Perform a PWM regression
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper})
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param motif Initial motif to start the regression from. The character "*" indicates a wildcard.
#' If NULL - a K-mer screen would be performed in order to find the best kmer for initialization.
#' @param motif_length Length of the seed motif. If the motif is shorter than this, it will be extended by wildcards (stars). Note that If the motif is longer than this, it will \emph{not} be truncated.
#' @param score_metric metric to use for optimizing the PWM. One of "r2" or "ks". For categorical response variables (0 and 1), "ks" is recommended, while for continuous response variables, "r2" is recommended. Default is "r2". When using "ks" the response variable should be a single vector of 0 and 1.
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
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{pssm: }{data frame with the pssm matrix with the inferred motif, where rows are positions and columns are nucleotides.}
#' \item{spat: }{a PSSM matrix with the inferred spatial model, where rows are positions and columns are nucleotides.}
#' \item{pred: }{a vector with the predicted pwm for each sequence.}
#' \item{response: }{The response matrix. If \code{include_response} is FALSE, the response matrix is not included in the list.}
#' \item{r2: }{\eqn{r^2} of the prediction with respect to the each response variable.}
#' \item{seed_motif: }{The seed motif that started the regression.}
#' \item(kmers: ){The k-mers that were screened in order to find the best seed motif (if motif was NULL).}
#' }
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#' res$pssm
#' res$spat
#' head(res$pred)
#'
#' plot_pssm_logo(res$pssm)
#'
#' @inheritParams screen_kmers
#' @inheritDotParams screen_kmers
#' @export
regress_pwm <- function(sequences,
                        response,
                        motif = NULL,
                        motif_length = 8,
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
            min_nuc_prob = min_nuc_prob,
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
        if (ncol(response) > 1 || (any(response[, 1] != 0 & response[, 1] != 1))) {
            cli_abort("When {.field score_metric} is {.val ks}, {.field response} should be a single vector of 0 and 1")
        }
    }

    cli_alert_info("Number of response variables: {.val {ncol(response)}}")


    if (is.null(motif)) {
        cli_alert_info("Screening for kmers in order to initialize regression")
        kmers <- screen_kmers(sequences, response, kmer_length = kmer_length, ...)
        motif <- kmers$kmer[which.max(abs(kmers$max_r2))]
        if (length(motif) == 0) { # could not find any kmer
            motif <- paste(rep("*", motif_length), collapse = "")
            cli_alert_info("Could not find any kmer. Initializing with {.val {motif}}")
        }
    } else {
        kmers <- NULL
    }

    if (stringr::str_length(motif) < motif_length) {
        motif <- stringr::str_pad(motif, motif_length, "*", side = "both")
        cli_alert_info("Motif is shorter than {.val {motif_length}}, extending with wildcards")
    }

    # replace gaps with wildcards
    stringr::str_replace(motif, "\\d+", function(x) {
        if (is.na(x)) {
            ""
        } else {
            paste(rep("*", x), collapse = "")
        }
    })

    cli_alert_info("Initializing regression with {.val {motif}}")

    cli_alert_info("Running regression")
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
        seed = seed
    )


    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2
    res$seed_motif <- motif

    if (!is.null(kmers)) {
        res$kmers <- kmers
    }

    cli_alert_success("Finished running regression")

    return(res)
}

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

#' Convert PSSM to kmer using majority
#'
#' @param pssm A PSSM matrix
#' @param min_freq minimal frequency of a nucleotide in the PSSM in order to be included in the kmer. Nuclotides with frequency less than this are set to "*".
#'
#' @return A kmer with the nucleotide with the highest frequency of each position in the PSSM. If there is no nucleotide with a high enough frequency, the nucleotide is set to "*".
#'
#'
#' @noRd
kmer_from_pssm <- function(pssm, min_freq = 0.3) {
    nucs <- c("A", "C", "G", "T")
    pssm <- pssm[, nucs]
    kmer <- apply(pssm, 1, function(x) {
        if (max(x) > min_freq) {
            return(nucs[which.max(x)])
        } else {
            return("*")
        }
    })
    kmer <- paste(kmer, collapse = "")
    return(kmer)
}
