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
#' @param max_gap The maximum length of a gap to be considered in the pattern. Default: 0
#' @param mask a string the length of \code{kmer_length} where 'N' indicates a wildcard position (default NULL - no mask).
#' @param add_mask if TRUE, the result of the mask will be added to the non-masked kmers. Otherwise - only the masked kmers would be returned.
#' @return A matrix where rows are the number of sequences, columns are the kmers and the values are the number of occurrences of each kmer.
#'
#' @examples
#' kmer_matrix(c("ATCG", "TCGA", "ATAT"), 2)
#' kmer_matrix(c("ATCG", "TCGA", "ATAT"), 3)
#' kmer_matrix(c("ATCG", "TCGA", "ATAT"), 3, mask = "ATN")
#'
#' @export
kmer_matrix <- function(sequences, kmer_length, max_gap = 0, mask = NULL, add_mask = FALSE, from_range = 1, to_range = NULL, set_rownames = FALSE) {
    if (length(sequences) == 0) {
        cli::cli_abort("{.field sequences} must have at least one element")
    }

    sequences <- toupper(sequences)

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

    if (!is.null(mask)) {
        if (stringr::str_length(mask) != kmer_length) {
            cli::cli_abort("{.field mask} must have the same length as {.field kmer_length}")
        }

        if (!any(grepl("N", mask))) {
            mask <- NULL
        }
    }

    mat <- kmer_matrix_cpp(sequences, kmer_length, from_range - 1, to_range, mask, add_mask, max_gap)
    if (set_rownames) {
        rownames(mat) <- sequences
    }
    return(mat)
}

#' Transform k-mers to PSSM (Position-Specific Scoring Matrix)
#'
#' This function transforms a vector of k-mers into a position-specific scoring matrix (PSSM).
#' A PSSM represents the frequency of each nucleotide at each position in the k-mers.
#' If a nucleotide is 'N', it is treated as equal probabilities for 'A', 'C', 'G', and 'T'.
#' The result is returned as a data frame with columns for the k-mer, position, and nucleotide frequencies.
#'
#' @param kmers A character vector of k-mers.
#' @param prior A numeric value indicating the prior probability for each nucleotide. Default is 0.01.
#'
#' @return A data frame with columns for the k-mer, position, and nucleotide frequencies, 'kmer', 'pos', 'A', 'C', 'G', 'T'.
#'
#' @examples
#' kmers_to_pssm(c("ACGTN", "TGCAN"), prior = 0.01)
#'
#' @export
kmers_to_pssm <- function(kmers, prior = 0.01) {
    if (length(grep("[^ACGTN]", kmers)) > 0) {
        cli::cli_abort("kmers must have only valid nucleotides")
    }

    calculate_pssm <- function(kmer) {
        # Initialize matrix
        pssm <- matrix(prior,
            nrow = 4, ncol = nchar(kmer),
            dimnames = list(c("A", "C", "G", "T"), 1:nchar(kmer))
        )

        # Fill matrix
        for (i in 1:nchar(kmer)) {
            if (substr(kmer, i, i) == "N") {
                pssm[, i] <- 1 / 4 # Equal probabilities for 'N'
            } else {
                pssm[substr(kmer, i, i), i] <- 1
            }
        }

        # renormalize
        pssm <- t(pssm)
        pssm <- pssm / rowSums(pssm)

        # Transform the matrix to data frame and reset row names
        pssm_df <- as.data.frame(pssm)
        colnames(pssm_df) <- c("A", "C", "G", "T")

        # Adding the 'pos' column
        pssm_df$pos <- 1:nchar(kmer)

        # Adding the 'kmer' column
        pssm_df$kmer <- kmer

        return(pssm_df)
    }

    # Apply function to each kmer and combine results
    pssm_df_all <- plyr::ldply(kmers, calculate_pssm, .id = "kmer") %>%
        select(kmer, pos, A, C, G, T)

    return(pssm_df_all)
}

#' Transform PSSM (Position-Specific Scoring Matrix) to a KMER
#'
#' This function transforms a PSSM into a k-mer of a given length.
#'
#' @param pssm PSSM matrix or data frame. The PSSM must have at least kmer_length rows.
#' @param kmer_length The length of the k-mer to return.
#' @param pos_bits_thresh A numeric value indicating the minimum number of bits per position to include the nucleotide in the k-mer. If the nucleotide does not meet this threshold, it is replaced with 'N'. Default is NULL.
#'
#' @return A character vector of length 1 containing the k-mer.
#'
#' @examples
#' pssm_to_kmer(get_motif_pssm("HOMER.AP_1"))
#' plot_pssm_logo_dataset("HOMER.AP_1")
#'
#' @export
pssm_to_kmer <- function(pssm, kmer_length = 7, pos_bits_thresh = NULL) {
    if (nrow(pssm) < kmer_length) {
        cli::cli_abort("pssm must have at least kmer_length rows")
    }

    bits <- bits_per_pos(pssm)
    bits[is.na(bits)] <- 0

    bits <- zoo::rollsum(bits, kmer_length, fill = NA, align = "left", na.rm = TRUE)
    pos <- which.max(bits)
    m <- pssm_to_mat(pssm)[pos:(pos + kmer_length - 1), ]
    kmer <- colnames(m)[apply(m, 1, which.max)]
    if (!is.null(pos_bits_thresh)) {
        bits <- bits_per_pos(m)
        kmer <- ifelse(bits > pos_bits_thresh, kmer, "N")
    }
    kmer <- paste(kmer, collapse = "")

    return(kmer)
}
