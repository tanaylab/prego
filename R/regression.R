#' Perform a PWM regression
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper})
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param motif Initial motif to start the regression from. The character "*" indicates a wildcard.
#' If NULL - a K-mer screen would be performed in order to find the best kmer for initialization.
#' @param motif_length Length of the seed motif. If the motif is shorter than this, it will be extended by wildcards (stars). Note that If the motif is longer than this, it will \emph{not} be truncated.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#' @param spat_min start of the spatial model from the beginning of the sequence (in bp)
#' @param spat_max end of the spatial model from the beginning of the sequence (in bp). If NULL - the spatial model
#' would end at the end of the sequence.
#' @param spat_bin size of the spatial bin (in bp).
#' @param min_nuc_prob minimum nucleotide probability in every iteration
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
                        bidirect = TRUE,
                        spat_min = 0,
                        spat_max = NULL,
                        spat_bin = 50,
                        min_nuc_prob = 0.001,
                        is_train = NULL,
                        include_response = TRUE,
                        seed = 60427,
                        verbose = FALSE,
                        kmer_length = 8,
                        ...) {
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

    cli_alert_info("Number of response variables: {.val {ncol(response)}}")


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
        is_bidirect = bidirect,
        verbose = verbose,
        seed = seed
    )

    if (include_response) {
        res$response <- response
    }

    res$r2 <- tgs_cor(response, as.matrix(res$pred))[, 1]^2
    res$seed_motif <- motif

    cli_alert_success("Finished running regression")


    return(res)
}

