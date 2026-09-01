#' Expected local PWM scores over a base frequency matrix
#'
#' @description
#' Score every motif in a database against a per-position base frequency matrix
#' at every start position. Where \code{\link{compute_local_pwm}} scores one
#' concrete sequence, this scores an ensemble of sequences summarised by its
#' per-position nucleotide distribution.
#'
#' @details
#' For a motif of length \eqn{L} placed at start \eqn{j}, with motif column
#' \eqn{p_l} and frequency column \eqn{q_{j+l}}, there are two ways to combine
#' each pair of distributions. They differ only in where the log sits:
#'
#' \describe{
#'   \item{\code{combine = "multiply"}}{\eqn{\sum_l \log(q_{j+l} \cdot p_l)}.
#'     Multiplies the per-column probabilities, i.e. the log of the expected
#'     likelihood. On a flat ensemble every motif gets the same floor,
#'     \eqn{L \log(0.25)}, so scores are comparable across motifs. This is
#'     exact only if the positions of the ensemble are independent - the usual
#'     PWM assumption, but one that an ensemble with covarying positions
#'     violates.}
#'   \item{\code{combine = "sum"}}{\eqn{\sum_l q_{j+l} \cdot \log p_l}. Sums the
#'     per-column log-probabilities, i.e. the expected log-likelihood - the mean
#'     score you would get by drawing sequences from the ensemble and running
#'     \code{\link{compute_local_pwm}} on each. Being linear in the frequencies,
#'     it is exact whatever the joint distribution of positions is. Its
#'     downside is that a flat ensemble gives each motif a different floor,
#'     strongly anti-correlated with the motif's information content, so scores
#'     are not comparable across motifs without normalisation.}
#' }
#'
#' Both reduce to \code{\link{compute_local_pwm}} when the frequency matrix is
#' one-hot, i.e. when the ensemble is a single sequence.
#'
#' @param freqs Base frequencies, in any of:
#'   \itemize{
#'     \item a single matrix, either \code{m x 4} or \code{4 x m}, in A,C,G,T
#'       order. Returns a single matrix rather than a list.
#'     \item a list of such matrices, which may differ in length.
#'     \item a 3-D array, \code{n_matrices x m x 4} or \code{n_matrices x 4 x m}.
#'     \item a stacked matrix, \code{n_matrices x 4m}, the layout
#'       \code{\link{calc_seq_pwm}} takes. To stack exactly 4 matrices use a
#'       list instead - a 4-row matrix is read as a single \code{4 x m} one.
#'   }
#'   Each position must be a distribution, i.e. sum to 1.
#' @param motifs A \code{MotifDB} object, or a tidy motif data frame with
#'   columns \code{motif}, \code{pos}, \code{A}, \code{C}, \code{G}, \code{T}.
#' @param combine How to combine the per-column products, \code{"multiply"} (the
#'   default) or \code{"sum"}. See Details.
#' @param bidirect Score both strands, combining them at each position as
#'   \code{log(exp(fwd) + exp(rev))} - the same convention as
#'   \code{\link{compute_local_pwm}}. \code{FALSE} scores the forward strand
#'   only.
#'
#' @return A numeric matrix with one row per motif and one column per position
#'   of the input, holding the score of the motif *starting* at that position.
#'   The last \code{L - 1} entries of each row are \code{NA}, where \code{L} is
#'   that motif's length. A list of such matrices if \code{freqs} held more than
#'   one.
#'
#' @examples
#' db <- all_motif_datasets()
#' mdb <- create_motif_db(db[db$motif %in% head(unique(db$motif), 20), ])
#'
#' # An ensemble that is certain at every position is just a sequence, and
#' # scores the same as compute_local_pwm() does on that sequence.
#' seq <- "ACGTACGTTGCAAGGTCCATACGT"
#' onehot <- diag(4)[match(strsplit(seq, "")[[1]], c("A", "C", "G", "T")), ]
#' scores <- calc_freq_local_pwm(onehot, mdb)
#' dim(scores)
#' scores[1:3, 1:5]
#'
#' # Under "multiply", a flat ensemble scores L * log(0.25) for every motif -
#' # the floor depends on the motif's length, not on the motif.
#' flat <- matrix(0.25, nrow = 24, ncol = 4)
#' calc_freq_local_pwm(flat, mdb, bidirect = FALSE)[1:3, 1]
#' mdb@motif_lengths[1:3] * log(0.25)
#'
#' # Several ensembles in one call
#' res <- calc_freq_local_pwm(list(certain = onehot, flat = flat), mdb)
#' sapply(res, dim)
#'
#' @seealso \code{\link{compute_local_pwm}}, \code{\link{calc_seq_pwm}}
#' @export
calc_freq_local_pwm <- function(freqs, motifs, combine = c("multiply", "sum"),
                                bidirect = TRUE) {
    combine <- match.arg(combine)

    # Keep the per-block BLAS dgemm serial; the TBB parallelFor over
    # (matrix, motif block) tasks is the parallelism layer.
    local_serial_blas()

    mdb <- as_motif_db(motifs)
    single <- is.matrix(freqs) && (ncol(freqs) == 4 || nrow(freqs) == 4)
    freq_list <- as_freq_list(freqs)

    res <- calc_freq_local_pwm_cpp(
        freq_list,
        mdb@mat,
        mdb@rc_mat,
        mdb@motif_lengths,
        multiply = combine == "multiply",
        bidirect = bidirect
    )

    res <- lapply(res, function(x) {
        rownames(x) <- colnames(mdb@mat)
        x
    })
    names(res) <- names(freq_list)

    if (single) res[[1]] else res
}

as_motif_db <- function(motifs) {
    if (methods::is(motifs, "MotifDB")) {
        return(motifs)
    }
    if (is.data.frame(motifs)) {
        return(create_motif_db(motifs))
    }
    cli_abort("{.field motifs} must be a {.cls MotifDB} object or a tidy motif data frame")
}

# Bring every accepted input form to a list of m x 4 matrices.
as_freq_list <- function(freqs) {
    if (is.list(freqs)) {
        if (length(freqs) == 0) {
            cli_abort("{.field freqs} is empty")
        }
        return(lapply(freqs, as_freq_matrix))
    }

    if (is.array(freqs) && length(dim(freqs)) == 3) {
        d <- dim(freqs)
        if (d[2] != 4 && d[3] != 4) {
            cli_abort("A 3-D {.field freqs} must be {.val n x m x 4} or {.val n x 4 x m}")
        }
        res <- lapply(seq_len(d[1]), function(k) as_freq_matrix(freqs[k, , ]))
        names(res) <- dimnames(freqs)[[1]]
        return(res)
    }

    if (!is.matrix(freqs)) {
        cli_abort("{.field freqs} must be a matrix, a list of matrices, or a 3-D array")
    }

    # A single matrix, in either orientation, takes precedence over the stacked
    # form; ncol == 4 first so that an m x 4 matrix with m divisible by 4 is
    # never read as stacked.
    if (ncol(freqs) == 4 || nrow(freqs) == 4) {
        return(list(as_freq_matrix(freqs)))
    }

    if (ncol(freqs) %% 4 != 0) {
        cli_abort("A stacked {.field freqs} must have a number of columns divisible by 4")
    }
    res <- lapply(seq_len(nrow(freqs)), function(k) {
        as_freq_matrix(matrix(freqs[k, ], ncol = 4, byrow = TRUE))
    })
    names(res) <- rownames(freqs)
    res
}

as_freq_matrix <- function(q) {
    q <- as.matrix(q)
    if (!is.numeric(q)) {
        cli_abort("{.field freqs} must be numeric")
    }
    if (ncol(q) != 4) {
        if (nrow(q) != 4) {
            cli_abort("Each frequency matrix must have 4 rows or 4 columns (A, C, G, T)")
        }
        q <- t(q)
    }

    sums <- rowSums(q)
    if (any(!is.finite(sums)) || any(abs(sums - 1) > 1e-6)) {
        hint <- NULL
        if (nrow(q) == 4 && ncol(q) %% 4 == 0) {
            hint <- c("i" = "A {.val 4 x 4m} matrix is read as a single {.val 4 x m} one. To pass 4 stacked matrices, use a list.")
        }
        cli_abort(c(
            "Every position of {.field freqs} must sum to 1.",
            "x" = "{sum(abs(sums - 1) > 1e-6)} position{?s} do{?es/} not.",
            hint
        ))
    }
    storage.mode(q) <- "double"
    unname(q)
}
