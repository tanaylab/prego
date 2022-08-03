
#' Screen for kmers
#'
#' @param sequences A vector of sequences (will go thruogh toupper)
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param kmer_length The number of non-gap characters in motifs that will be screened
#' @param min_cor Only patterns for which the maximum correlation to one of the response variable is larger than min_cor will be reported
#' @param min_n Only patterns for which the average number of occurrences in the sequences is larger than min_n will be reported
#' @param is_train a boolean vector that determine which subset of sequences to use when screening
#' @param min_gap,max_gap the length of a gap to be considered in the pattern. Only one gap, of length min_gap:max_gap, is being used, and is located anywhere in the motif. Note that this greatly expand the search space (and increase multiple tesing severly).
#' @param from_range Sequences will be considered only from position from_range (default 0)
#' @param to_range Sequences will be considered only up to position to_range (default NULL - using the length of the sequences)
#' @param seed random seed
#'
#' @return A data frame with the following columns, together with a column for each response variable with the
#' correlation of the kmers to the response variable:
#' \itemize{
#' \item{kmer:}{the kmer pattern},
#' \item{max_r2:}{the maximum R^2 to one of the response variables},
#' \item{avg_n:}{the average number of times the kmer appears in the sequences},
#' \item{avg_var:}{the variance of the number of times the kmer appears in the sequences},
#' }
#'
#'
#'
#'
#' @export
screen_kmers <- function(sequences,
                         response,
                         kmer_length = 6,
                         min_cor = 0.08,
                         min_n = 50,
                         is_train = NULL,
                         min_gap = 0,
                         max_gap = 0,
                         from_range = 0,
                         to_range = NULL,
                         seed = 60427) {
    if (is.null(is_train)) {
        is_train <- rep(TRUE, length(sequences))
    }

    n_in_train <- sum(is_train)

    if (is.null(to_range)) {
        to_range <- nchar(sequences[1])
    }

    res <- screen_kmers_cpp(
        toupper(sequences),
        response,
        is_train,
        kmer_length,
        from_range,
        to_range,
        min_cor,
        min_n,
        min_gap,
        max_gap,
        n_in_train,
        seed
    )

    return(res)
}
