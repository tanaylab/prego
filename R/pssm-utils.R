#' Compute PWMs for a set of sequences given a PSSM matrix
#'
#' @param sequences a vector of sequences
#' @param pssm a PSSM matrix or data frame. The columns of the matrix or data frame should be named with the nucleotides ('A', 'C', 'G' and 'T').
#' @param spat a data frame with the spatial model (as returned from the \code{$spat} slot from the regression). Should contain a column called 'bin' and a column called 'spat_factor'.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#'
#' @return a vector with the predicted pwm for each sequence.
#'
#' @examples
#' res <- regress_pwm(sequences_example, response_mat_example)
#'
#' pwm <- compute_pwm(sequences_example, res$pssm)
#' head(pwm)
#'
#' # this is similar to the prediction in the regression
#' head(res$pred)
#'
#' @inheritParams regress_pwm
#' @export
compute_pwm <- function(sequences, pssm, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE) {
    if (is.null(spat)) {
        spat <- data.frame(bin = 0, spat_factor = 1)
        binsize <- nchar(sequences[[1]])
    } else {
        validate_spat(spat)
        binsize <- unique(diff(spat$bin))
    }

    if (is.null(spat_max)) {
        spat_max <- nchar(sequences[1])
    }

    if (!all(c("A", "C", "G", "T") %in% colnames(pssm))) {
        cli_abort("The {.field pssm} matrix should have columns {.val A}, {.val C}, {.val G}, {.val T}")
    }

    pwm <- compute_pwm_cpp(
        sequences = toupper(sequences),
        pssm = as.matrix(pssm[, c("A", "C", "G", "T")]),
        is_bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_factor = spat$spat_factor,
        bin_size = binsize
    )

    return(pwm)
}

validate_spat <- function(spat) {
    if (!is.data.frame(spat)) {
        cli_abort("The {.field spat} argument should be a data frame")
    }
    if (!all(c("bin", "spat_factor") %in% colnames(spat))) {
        cli_abort("The {.field spat} data frame should have columns {.val bin} and {.val spat_factor}")
    }
    if (!is.numeric(spat$bin) || !is.numeric(spat$spat_factor)) {
        cli_abort("The {.field spat} data frame should have columns {.val bin} and {.val spat_factor} of type numeric")
    }
    binsize <- unique(diff(spat$bin))
    if (length(binsize) > 1) {
        cli_abort("The bins in {.field spat} should be of equal size")
    }
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
