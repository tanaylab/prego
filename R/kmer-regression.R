#' Generate kmers
#'
#' This function generates all possible kmers considering the gap length. Gaps are represented by 'N'.
#'
#' @param kmer_length The number of non-gap characters in motifs that will be screened.
#' @param max_gap The maximum length of a gap to be considered in the pattern. Default: 0
#' @param min_gap The minimum length of a gap to be considered in the pattern. Default: 0
#' @return A vector of all possible kmers considering the gap length.
#'
#' @examples
#'
#' # Generate kmers of length 2 without any gaps
#' generate_kmers(2)
#'
#' # Generate kmers of length 3 with a single gap (1 'N') at any position
#' generate_kmers(3, min_gap = 1, max_gap = 1)
#'
#' # Generate kmers of length 3 with a gap of 1 to 2 'N's at any position
#' generate_kmers(3, min_gap = 1, max_gap = 2)
#'
#' # Generate kmers of length 3 with a gap of 2 'N's at any position
#' generate_kmers(3, min_gap = 2, max_gap = 2)
#'
#' @export
generate_kmers <- function(kmer_length, max_gap = 0, min_gap = 0) {
    dna <- c("T", "C", "G", "A")

    if (kmer_length < 1) {
        cli::cli_abort("{.field kmer_length} must be greater than 0")
    }

    if (min_gap < 0) {
        cli::cli_abort("{.field min_gap} must be greater than or equal to 0")
    }

    if (max_gap < 0) {
        cli::cli_abort("{.field max_gap} must be greater than or equal to 0")
    }

    if (max_gap > kmer_length) {
        cli::cli_abort("{.field max_gap} must be less than or equal to {.field kmer_length}")
    }

    if (max_gap < min_gap) {
        cli::cli_abort("{.field max_gap} must be greater than or equal to {.field min_gap}")
    }

    gap <- seq(min_gap, max_gap)
    kmers <- expand.grid(rep(list(dna), kmer_length))
    kmers <- apply(kmers, 1, paste, collapse = "")

    gap_kmers <- c()

    for (g in gap) {
        for (pos in 1:(kmer_length - g + 1)) {
            gap_seq <- paste(rep("N", g), collapse = "")
            gap_kmers <- c(gap_kmers, paste(substr(kmers, 1, pos - 1), gap_seq, substr(kmers, pos + g, kmer_length), sep = ""))
        }
    }

    if (min_gap > 0) {
        return(gap_kmers)
    }

    return(unique(c(kmers, gap_kmers)))
}

#' Generate a kmer Matrix
#'
#' This function calculates the frequency of each kmer for each DNA sequence.
#'
#' @param sequences A vector of strings with DNA sequences ('T', 'C', 'G', 'A' or 'N').
#' @param from_range Sequences will be considered only from position from_range.
#' @param to_range Sequences will be considered only up to position to_range (default NULL - using the length of the sequences).
#' @param set_rownames If TRUE, the rownames of the matrix will be set to the sequences (default FALSE).
#' @return A matrix where rows are the number of sequences, columns are the possible kmers and the values are the number of occurrences of each kmer.
#'
#' @examples
#' kmer_matrix(c("ATCG", "ATCG"), 2)
#' kmer_matrix(c("ATCG", "ATCG"), 2, 1)
#'
#' @inheritParams generate_kmers
#' @export
kmer_matrix <- function(sequences, kmer_length, max_gap = 0, min_gap = 0, from_range = 1, to_range = NULL, set_rownames = FALSE) {
    if (length(sequences) == 0) {
        cli::cli_abort("{.field sequences} must have at least one element")
    }

    if (from_range < 1) {
        cli::cli_abort("{.field from_range} must be greater than 0")
    }

    if (!is.null(to_range)) {
        if (to_range < from_range) {
            cli::cli_abort("{.field to_range} must be greater than or equal to {.field from_range}")
        }

        if (to_range > nchar(sequences[1])) {
            cli::cli_abort("{.field to_range} must be less than or equal to the length of the sequences")
        }
    }


    kmers <- generate_kmers(kmer_length = kmer_length, max_gap = max_gap, min_gap = min_gap)
    mat <- kmer_matrix_cpp(sequences, kmers, from_range - 1, to_range)
    if (set_rownames) {
        rownames(mat) <- sequences
    }
    return(mat)
}
