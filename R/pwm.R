num_to_one_hot <- function(x, bits = 4L) {
    as.vector(diag(1L, bits)[, x])
}

seqs_to_onehot <- function(seqs) {
    # Convert sequences to one-hot encoded matrix
    if (is.vector(seqs)) {
        seqs <- as.list(seqs) # Handle single sequence passed as vector
    } else if (is.character(seqs)) {
        seqs <- list(seqs) # Handle single sequence passed as string
    }

    # Split and convert all sequences at once
    onehot_list <- lapply(seqs, function(s) {
        strsplit(s, "")[[1]] |>
            match(c("A", "C", "G", "T")) |>
            num_to_one_hot(bits = 4L)
    })

    # Combine into single matrix
    onehot_mat <- do.call(rbind, onehot_list)

    return(onehot_mat)
}

#' Calculate Position Weight Matrix (PWM) Scores for DNA Sequences
#'
#' @param sequences Character vector of DNA sequences.
#' @param mdb MotifDB object containing PWMs.
#'
#' @return A numeric matrix with sequences as rows and motifs as columns, containing PWM scores.
#'   Row names are preserved from input sequences if they exist.
#'   Column names are preserved from the PWM matrix if they exist.
#'
#' @inheritParams compute_pwm
#' @examples
#' sequences <- c("ACGTACGT", "TGCATGCA")
#' pwm <- matrix(runif(16), nrow = 4)
#' scores <- calc_seq_pwm(sequences, pwm)
#'
#' @export
calc_seq_pwm <- function(sequences, mdb, bidirect = TRUE) {
    # Input validation
    if (!is.character(sequences)) {
        stop("sequences must be a character vector")
    }

    # Convert sequences to uppercase and save original names
    sequences <- toupper(sequences)
    seq_names <- names(sequences)
    pwm_names <- colnames(mdb@mat)

    # Convert sequences to one-hot encoding
    onehot_seqs <- seqs_to_onehot(sequences)

    # Calculate PWM scores with spatial factors
    result <- calc_seq_pwm_parallel_cpp(
        onehot_seqs,
        mdb@mat,
        mdb@rc_mat,
        mdb@motif_lengths,
        min(mdb@motif_lengths),
        bidirect,
        mdb@spat_factors,
        mdb@spat_bin_size
    )

    # Restore names if they existed
    if (!is.null(seq_names)) {
        rownames(result) <- seq_names
    }
    if (!is.null(pwm_names)) {
        colnames(result) <- pwm_names
    }

    return(result)
}
