#' Perform a PWM regression
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper})
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param motif Initial motif to start the regression from. If NULL - a K-mer screen would be performed in order
#' to find the best kmer for initialization.
#' @param verbose show verbose messages.
#' @param seed random seed
#'
#' @inheritParams screen_kmers
#' @inheritDotParams screen_kmers
#' @export
regress_pwm <- function(sequences,
                        response,
                        motif = NULL,
                        bidirect = TRUE,
                        epsilon = 0.001,
                        spat_min = 0,
                        spat_max = NULL,
                        min_nuc_prob = 0.001,
                        spat_bin = 50,
                        min_rms_for_star = 0.001,
                        is_train = NULL,
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
            motif <- paste(rep("*", kmer_length), collapse = "")
            cli_alert_info("Could not find any kmer. Initializing with {.val {motif}}")
        }
    }
    cli_alert_info("Initializing regression with {.val {motif}}")
    cli_alert_info("Running regression")

    res <- regress_pwm_cpp(
        toupper(sequences),
        response,
        is_train,
        motif = motif,
        epsilon = epsilon,
        min_rms_for_star = min_rms_for_star,
        spat_min = spat_min,
        spat_max = spat_max,
        min_nuc_prob = min_nuc_prob,
        spat_bin = spat_bin,
        is_bidirect = bidirect,
        verbose = verbose,
        seed = seed
    )

    cli_alert_success("Finished running regression")

    return(res)
}

#' Plot LOGO of the pssm result from the regression
#'
#' @param pssm the 'pssm' field from the regression result
#'
#'
#' @export
plot_pssm_logo <- function(pssm) {
    pfm <- pssm %>%
        select(-1) %>%
        as.data.frame() %>%
        column_to_rownames("pos") %>%
        as.matrix() %>%
        t()
    ggseqlogo::ggseqlogo(pfm)
}
