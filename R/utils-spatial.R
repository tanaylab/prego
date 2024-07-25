calculate_bins <- function(max_seq_len, spat_num_bins = NULL, spat_bin_size = NULL, default_bin_size = 40) {
    if (!is.null(spat_num_bins) && !is.null(spat_bin_size)) {
        return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
    }

    if (!is.null(spat_num_bins)) {
        spat_bin_size <- floor(max_seq_len / spat_num_bins)
        return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
    }

    if (!is.null(spat_bin_size)) {
        spat_num_bins <- floor(max_seq_len / spat_bin_size)
        if (spat_num_bins %% 2 == 0) {
            spat_num_bins <- spat_num_bins - 1
        }
        return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
    }

    spat_bin_size <- default_bin_size
    spat_num_bins <- floor(max_seq_len / spat_bin_size)

    if (spat_num_bins %% 2 == 0) {
        spat_num_bins <- spat_num_bins - 1
    }

    if (spat_num_bins < 3) {
        cli_abort("Calculated spat_num_bins is less than 3, which is not permissible.")
    }

    return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
}

calc_spat_min_max <- function(spat_num_bins, max_seq_len, spat_bin_size) {
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

#' Calculate Dinucleotide Distribution in Sequences
#'
#'
#' @param sequences a character vector containing the sequences to analyze.
#'        Each element of the vector should be a single sequence.
#' @param size an integer specifying the size to consider for the analysis.
#'        If NULL (default), the maximum length of the sequences in the `sequences`
#'        vector is used.
#'
#' @return a data frame with columns 'pos' and 16 columns representing each possible
#'         dinucleotide. Each row represents a position in the sequences (from 1 to `size`),
#'         and contains the fraction of each dinucleotide at that position across
#'         all sequences.
#'
#' @examples
#'
#' # Generate some random sequences for testing
#' set.seed(60427)
#' sequences <- sapply(1:100, function(x) {
#'     paste0(sample(c("A", "C", "G", "T"), 1000, replace = TRUE), collapse = "")
#' })
#' sequences <- as.character(sequences)
#'
#' # Calculate the dinucleotide distribution
#' result <- calc_sequences_dinuc_dist(sequences)
#'
#' head(result)
#'
#' @export
calc_sequences_dinuc_dist <- function(sequences, size = NULL) {
    return(calc_sequences_n_nuc_dist(sequences, n = 2, size = size))
}

#' Calculate Dinucleotide Counts for Sequences
#'
#' This function calculates the total count of each dinucleotide for each sequence
#' in a vector of DNA sequences.
#'
#' @param sequences A character vector of DNA sequences. Each element should be
#'   a string representing a DNA sequence composed of A, T, C, and G.
#'
#' @return A numeric matrix where:
#'   * Each row corresponds to a sequence in the input vector.
#'   * Each column represents a specific dinucleotide (AA, AC, AG, AT, CA, CC, etc.).
#'   * The values in the matrix are the counts of each dinucleotide in each sequence.
#'   * Column names are set to the corresponding dinucleotides.
#'
#' @examples
#'
#' sequences <- c("ATCG", "GCTA", "AATT")
#' result <- calc_sequences_dinucs(sequences)
#' print(result)
#'
#' @export
calc_sequences_dinucs <- function(sequences) {
    res <- calc_sequences_dinuc_cpp(toupper(sequences))
    if (!is.null(names(sequences))) {
        rownames(res) <- names(sequences)
    }
    return(res)
}

#' Calculate Trinucleotide Distribution in Sequences
#'
#' @param sequences a character vector containing the sequences to analyze.
#'        Each element of the vector should be a single sequence.
#' @param size an integer specifying the size to consider for the analysis.
#'        If NULL (default), the maximum length of the sequences in the `sequences`
#'        vector is used.
#'
#' @return a data frame with columns 'pos' and 64 columns representing each possible
#'         trinucleotide. Each row represents a position in the sequences (from 1 to `size`),
#'         and contains the fraction of each trinucleotide at that position across
#'         all sequences.
#'
#' @examples
#'
#' # Generate some random sequences for testing
#' set.seed(60427)
#' sequences <- sapply(1:100, function(x) {
#'     paste0(sample(c("A", "C", "G", "T"), 1000, replace = TRUE), collapse = "")
#' })
#' sequences <- as.character(sequences)
#'
#' # Calculate the trinucleotide distribution
#' result <- calc_sequences_trinuc_dist(sequences)
#'
#' head(result)
#'
#' @export
calc_sequences_trinuc_dist <- function(sequences, size = NULL) {
    return(calc_sequences_n_nuc_dist(sequences, n = 3, size = size))
}

#' Calculate n-nucleotide Distribution in Sequences
#'
#' @param sequences a character vector containing the sequences to analyze.
#'        Each element of the vector should be a single sequence.
#' @param n an integer specifying the length of the n-nucleotide (e.g., 2 for dinucleotides, 3 for trinucleotides).
#' @param size an integer specifying the size to consider for the analysis.
#'        If NULL (default), the maximum length of the sequences in the `sequences`
#'        vector is used.
#'
#' @return a data frame with columns 'pos' and 4^n columns representing each possible
#'         n-nucleotide. Each row represents a position in the sequences (from 1 to `size`),
#'         and contains the fraction of each n-nucleotide at that position across
#'         all sequences.
#'
#' @examples
#'
#' # Generate some random sequences for testing
#' set.seed(60427)
#' sequences <- sapply(1:100, function(x) {
#'     paste0(sample(c("A", "C", "G", "T"), 1000, replace = TRUE), collapse = "")
#' })
#' sequences <- as.character(sequences)
#'
#' # Calculate the n-nucleotide distribution (e.g., for trinucleotides)
#' result <- calc_sequences_n_nuc_dist(sequences, n = 3)
#'
#' head(result)
#'
#' @noRd
calc_sequences_n_nuc_dist <- function(sequences, n, size = NULL) {
    if (!is.character(sequences)) {
        cli::cli_abort("The sequences must be a character vector")
    }

    if (!is.numeric(n) || n < 1) {
        cli::cli_abort("n must be a positive integer")
    }

    if (is.null(size)) {
        size <- max(nchar(sequences))
    }

    if (any(nchar(sequences) < size)) {
        cli::cli_abort("Some sequences are shorter than the size")
    }

    result <- n_nuc_distribution(sequences, n = n, size = size)

    # set last n-1 positions to NA
    result[(size - n + 2):size, -1] <- NA

    return(result)
}
