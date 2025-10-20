#' Compute PWMs for a set of sequences given a PSSM matrix
#'
#' @param sequences a vector of sequences
#' @param pssm a PSSM matrix or data frame. The columns of the matrix or data frame should be named with the nucleotides ('A', 'C', 'G' and 'T').
#' @param spat a data frame with the spatial model (as returned from the \code{$spat} slot from the regression). Should contain a column called 'bin' and a column called 'spat_factor'.
#' @param spat_min the minimum position to use from the sequences. The default is 1.
#' @param spat_max the maximum position to use from the sequences. The default is the length of the sequences.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#' @param prior a prior probability for each nucleotide.
#' @param func the function to use to combine the PWMs for each sequence. Either 'logSumExp' or 'max'. The default is 'logSumExp'.
#'
#' @return a vector with the predicted pwm for each sequence.
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#'
#' pwm <- compute_pwm(cluster_sequences_example, res$pssm, res$spat)
#' head(pwm)
#'
#' # this is similar to the prediction in the regression
#' head(res$pred)
#' }
#'
#' @export
compute_pwm <- function(sequences, pssm, spat = NULL, spat_min = 1, spat_max = NULL, bidirect = TRUE, prior = 0.01, func = "logSumExp") {
    if (is.null(spat)) {
        spat <- data.frame(bin = 0, spat_factor = 1)
        binsize <- nchar(sequences[[1]])
    } else {
        validate_spat(spat)
        binsize <- unique(diff(spat$bin))
    }

    if (is.null(spat_max) || is.na(spat_max)) {
        spat_max <- nchar(sequences[[1]])
    }

    if (is.null(spat_min) || is.na(spat_max)) {
        spat_min <- 1
    }

    if (!(spat_min == 1 && spat_max == nchar(sequences[[1]]))) {
        sequences <- stringr::str_sub(sequences, start = spat_min, end = spat_max)
    }

    if (!all(c("A", "C", "G", "T") %in% colnames(pssm))) {
        cli_abort("The {.field pssm} matrix should have columns {.val A}, {.val C}, {.val G}, {.val T}")
    }

    pssm_mat <- as.matrix(pssm[, c("A", "C", "G", "T")])
    pssm_mat <- pssm_mat / rowSums(pssm_mat)

    if (prior < 0 || prior > 1) {
        cli_abort("The {.field prior} should be between 0 and 1")
    }

    if (prior > 0) {
        pssm_mat <- pssm_mat + prior
        pssm_mat <- pssm_mat / rowSums(pssm_mat)
    }

    if (func == "max") {
        use_max <- TRUE
    } else if (func == "logSumExp") {
        use_max <- FALSE
    } else {
        cli_abort("The {.field func} argument should be either {.val max} or {.val logSumExp}")
    }

    pwm <- compute_pwm_cpp(
        sequences = toupper(sequences),
        pssm_mat = pssm_mat,
        is_bidirect = bidirect,
        spat_min = 0,
        spat_max = nchar(sequences[1]),
        spat_factor = spat$spat_factor,
        bin_size = binsize,
        use_max = use_max
    )

    return(pwm)
}

#' Compute local PWMs for a set of sequences given a PSSM matrix
#'
#' @description compute the local PWM for each position in every sequence. The edges of each sequences would become NA.
#'
#' @param return_list Logical. If TRUE, returns a list with one vector per sequence (useful for sequences of different lengths). If FALSE, returns a matrix (requires all sequences to have the same length). Default is FALSE.
#'
#' @return If return_list is FALSE: a matrix with \code{length(sequences)} rows and \code{max(nchar(sequences))} columns with the local PWM for each sequence in each position. If return_list is TRUE: a list with \code{length(sequences)} elements, where each element is a numeric vector of local PWM scores for the corresponding sequence.
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#'
#' # Return matrix (sequences must have same length)
#' pwm <- compute_local_pwm(cluster_sequences_example, res$pssm, res$spat)
#' head(pwm)
#'
#' # Return list (allows sequences of different lengths)
#' pwm_list <- compute_local_pwm(cluster_sequences_example, res$pssm, res$spat, return_list = TRUE)
#' pwm_list[[1]]
#'
#' # Using a motif from MOTIF_DB
#' hnf1a_pwm <- compute_local_pwm(cluster_sequences_example, as.matrix(MOTIF_DB["JASPAR.HNF1A"]), return_list = TRUE)
#' hnf1a_pwm[[1]]
#' }
#'
#' @inheritParams compute_pwm
#' @export
compute_local_pwm <- function(sequences, pssm, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01, return_list = FALSE) {
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

    pssm_mat <- as.matrix(pssm[, c("A", "C", "G", "T")])

    if (nrow(pssm_mat) == 0) {
        cli_abort("The {.field pssm} matrix should have at least one row")
    }

    if (prior < 0 || prior > 1) {
        cli_abort("The {.field prior} should be between 0 and 1")
    }

    if (prior > 0) {
        pssm_mat <- pssm_mat + prior
    }

    # Check sequence lengths only if not returning a list
    if (!return_list) {
        seq_l <- stringr::str_length(sequences)
        if (length(unique(seq_l)) > 1) {
            cli_abort("All sequences should have the same length when {.field return_list} is FALSE. Set {.field return_list = TRUE} to handle sequences of different lengths.")
        }
    }

    pwm <- compute_local_pwm_cpp(
        sequences = toupper(sequences),
        pssm_mat = pssm_mat,
        is_bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_factor = spat$spat_factor,
        bin_size = binsize,
        return_list = return_list
    )

    return(pwm)
}

#' Calculate Theoretical Scores for a Position-Specific Scoring Matrix (PSSM)
#'
#' These functions compute theoretical scores (maximum, minimum, and quantiles) that can be achieved
#' by sequences matching a given Position-Specific Scoring Matrix (PSSM).
#'
#' @param pssm A matrix or data frame containing position-specific probabilities for nucleotides A, C, G, T
#' @param prior A numeric value (default: 0.01) added to each probability to avoid zero probabilities
#' @param regularization A numeric value (default: 0.01) added inside the log function to prevent -Inf values
#' @param q A numeric value between 0 and 1 specifying the quantile to calculate (only for pssm_quantile)
#'
#' @return
#' \itemize{
#'   \item pssm_theoretical_max: Maximum possible score for the given PSSM
#'   \item pssm_theoretical_min: Minimum possible score for the given PSSM
#'   \item pssm_quantile: Score at the specified quantile between min and max scores
#' }
#'
#' @details The functions normalize the probabilities and calculate scores based on logarithms
#' with regularization. The quantile function interpolates linearly between min and max scores.
#'
#' @examples
#' pssm <- as.matrix(MOTIF_DB["HOMER.CTCF"])
#' pssm_theoretical_max(pssm)
#' pssm_theoretical_min(pssm)
#' pssm_quantile(pssm, 0.85)
#'
#' @export
pssm_theoretical_max <- function(pssm, prior = 0.01, regularization = 0.01) {
    pssm <- as.matrix(pssm[, c("A", "C", "G", "T")])
    if (prior > 0) {
        pssm <- pssm + prior
    }
    pssm <- pssm / rowSums(pssm)
    theo_max <- sum(log(regularization + apply(pssm, 1, max)))
    return(theo_max)
}

#' @rdname pssm_theoretical_max
#' @export
pssm_theoretical_min <- function(pssm, prior = 0.01, regularization = 0.01) {
    pssm <- as.matrix(pssm[, c("A", "C", "G", "T")])
    if (prior > 0) {
        pssm <- pssm + prior
    }
    pssm <- pssm / rowSums(pssm)
    theo_min <- sum(log(regularization + apply(pssm, 1, min)))
    return(theo_min)
}

#' @rdname pssm_theoretical_max
#' @export
pssm_quantile <- function(pssm, q, prior = 0.01, regularization = 0.01) {
    theo_max <- pssm_theoretical_max(pssm, prior, regularization)
    theo_min <- pssm_theoretical_min(pssm, prior, regularization)
    return(theo_min + q * (theo_max - theo_min))
}

#' Extracts local position weight matrix (PWM) scores for given intervals and a PWM.
#'
#' @param intervals The intervals to extract
#'
#' @return A matrix with \code{nrow(intervals)} rows and \code{ncol(pssm)} columns with the local PWM for each sequence in each position.
#'
#' @inheritParams compute_local_pwm
#' @export
gextract.local_pwm <- function(intervals, pssm, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01) {
    sequences <- intervals_to_seq(intervals)

    res <- compute_local_pwm(sequences = sequences, pssm = pssm, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior)
    rnames <- paste0(intervals$chrom, "_", intervals$start, "_", intervals$end)
    if (length(unique(rnames)) != length(rnames)) {
        rnames <- paste0(rnames, ".", 1:nrow(res))
    }
    rownames(res) <- rnames

    return(res)
}

#' Center intervals by PSSM
#'
#' This function takes a set of intervals and a position-specific scoring matrix (PSSM) and centers the intervals
#' based on the maximum score position in the PSSM. The intervals are shifted so that the maximum score position
#' becomes the center of each interval.
#'
#' @inheritParams gextract.local_pwm
#'
#' @return A data frame containing the centered intervals. The intervals will have the same columns as the input
#'   intervals, but the start and end positions will be adjusted to center the intervals based on the maximum score
#'   position in the PSSM.
#'
#' @export
gintervals.center_by_pssm <- function(intervals, pssm, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01) {
    local_pwm <- gextract.local_pwm(intervals, pssm, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior)
    maxs <- apply(local_pwm, 1, which.max)
    intervals_size <- intervals$end[1] - intervals$start[1]
    intervs_center <- intervals %>%
        mutate(start = start + maxs, end = start + 1) %>%
        misha.ext::gintervals.normalize(intervals_size) %>%
        select(chrom, start, end, everything())
    return(intervs_center)
}

#' Calculate the frequency of a position weight matrix (PWM) in a given set of intervals
#'
#' @param intervals The intervals to extract
#' @param q_threshold The quantile threshold of the PWM (e.g. 0.99 for the top percentile)
#' @inheritParams gextract.local_pwm
#' @inheritParams gpwm_quantiles
#'
#' @return a matrix with \code{nrow(intervals)} rows and \code{ncol(pssm)} columns with the TRUE if the PWM is above the threshold for each sequence in each position.
#'
#' @export
gextract.local_pwm_freq <- function(intervals, pssm, q_threshold, bg_intervals = NULL, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01, n_sequences = 1e4, dist_from_edge = 3e6, chromosomes = NULL) {
    local_pwm <- gextract.local_pwm(intervals, pssm, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior)
    size <- intervals$end[1] - intervals$start[1]
    if (is.null(bg_intervals)) {
        bg_intervals <- misha.ext::grandom_genome(size, n_sequences, dist_from_edge, chromosomes) %>%
            distinct(chrom, start, end)
    }
    local_pwm_r <- gextract.local_pwm(bg_intervals, pssm, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior)

    thresh <- quantile(local_pwm_r, q_threshold, na.rm = TRUE)
    pwm_freq <- local_pwm > thresh
    pwm_freq <- pwm_freq + 0

    return(pwm_freq)
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

#' Calculate the number of bits per position in a Position-Specific Scoring Matrix (PSSM).
#'
#' This function takes a PSSM as input and calculates the number of bits per position.
#' The PSSM should be a data frame or matrix with columns representing the nucleotides A, C, G, and T.
#' The function first normalizes the PSSM by dividing each element by the sum of its row.
#' Then, it calculates the entropy for each position using the formula: bits = log2(4) + sum(p * log2(p)),
#' where p is the probability of each nucleotide at the position.
#' Finally, it sets any negative values to zero and returns the resulting bits per position.
#'
#' @param pssm A data frame or matrix representing the Position-Specific Scoring Matrix (PSSM).
#' @param prior A numeric value indicating the prior probability for each nucleotide. Default is 0.01.
#' @return A numeric vector representing the number of bits per position in the PSSM.
#' @examples
#' pssm <- data.frame(
#'     A = c(0.2, 0.3, 0.1, 0.4),
#'     C = c(0.1, 0.2, 0.3, 0.4),
#'     G = c(0.4, 0.3, 0.2, 0.1),
#'     T = c(0.3, 0.2, 0.4, 0.1)
#' )
#' bits_per_pos(pssm)
#'
#' @export
bits_per_pos <- function(pssm, prior = 0.01) {
    pssm <- as.matrix(pssm[, c("A", "C", "G", "T")])
    if (prior > 0) {
        pssm <- pssm + prior
    }
    pssm <- pssm / rowSums(pssm)
    bits <- log2(4) + rowSums(pssm * log2(pssm))
    bits <- pmax(bits, 0)
    return(bits)
}

#' Mask sequences by thresholding the PWM
#'
#' @description Mask sequences by thresholding the PWM. Sequences with a PWM above the threshold will be masked by 'N'.
#' Sequences at the edges of the sequences will also be masked by 'N'.
#'
#' @param mask_thresh Threshold for masking. Sequences with a PWM above this threshold will be masked by 'N'.
#' @param pos_bits_thresh Mask only positions with amount of information contributed (Shannon entropy, measured in bits) above this threshold.
#' The scale is the same as the y axis in the pssm logo plots.
#'
#' @return A vector with the masked sequences.
#'
#' @examples
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' new_sequences <- mask_sequences_by_pwm(
#'     cluster_sequences_example,
#'     res$pssm,
#'     quantile(res$pred, 0.95),
#'     spat = res$spat
#' )
#'
#' head(new_sequences)
#'
#' @inheritParams compute_pwm
#' @export
mask_sequences_by_pwm <- function(sequences, pssm, mask_thresh, pos_bits_thresh = 0.2, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01) {
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

    pssm_mat <- as.matrix(pssm[, c("A", "C", "G", "T")])

    if (prior < 0 || prior > 1) {
        cli_abort("The {.field prior} should be between 0 and 1")
    }

    if (prior > 0) {
        pssm_mat <- pssm_mat + prior
    }

    bits <- bits_per_pos(pssm)
    pos_mask <- bits > pos_bits_thresh
    if (sum(pos_mask) == 0) {
        cli_warn("No positions with information content above {.val {pos_bits_thresh}} were found")
    } else {
        cli_alert_info("The following positions will be masked: {.val {which(pos_mask)}}. Overall {.val {sum(pos_mask)}} positions will be masked")
    }

    res <- mask_sequences_cpp(
        sequences = toupper(sequences),
        pssm_mat = pssm_mat,
        is_bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_factor = spat$spat_factor,
        bin_size = binsize,
        mask_thresh = mask_thresh,
        pos_mask = pos_mask
    )

    return(res)
}


#' Convert PSSM to consensus sequence
#'
#' @param pssm A PSSM matrix
#' @param single_thresh,double_thresh thresholds for the consensus sequence calculation
#' (single and double nucleotides)
#'
#' @return A consensus sequence for the PSSM. If no consensus sequence can be found, the function returns NA.
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' consensus_from_pssm(res$pssm)
#' }
#'
#' @export
consensus_from_pssm <- function(pssm, single_thresh = 0.4, double_thresh = 0.6) {
    consensus <- get_consensus_cpp(as.matrix(pssm[, c("A", "C", "G", "T")]), single_thresh, double_thresh)
    consensus <- gsub("^\\*+", "", gsub("\\*+$", "", consensus))
    if (consensus == "") {
        consensus <- NA
    }
    return(consensus)
}


pssm_to_mat <- function(pssm_df) {
    if (is.matrix(pssm_df)) {
        pssm_df <- pssm_df[, c("A", "C", "G", "T")]
        return(pssm_df)
    }

    if (!is.null(rownames(pssm_df))) {
        rownames(pssm_df) <- NULL
    }

    if (!("pos" %in% colnames(pssm_df))) {
        pssm_df$pos <- 1:nrow(pssm_df)
    }

    pssm_df %>%
        arrange(as.numeric(pos)) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("pos") %>%
        select(A, C, G, T) %>%
        as.matrix()
}

mat_to_pssm <- function(pssm) {
    as.data.frame(pssm) %>%
        mutate(pos = 1:nrow(pssm)) %>%
        select(pos, A, C, G, T)
}

pssm_add_prior <- function(pssm_df, prior) {
    pssm_mat <- pssm_to_mat(pssm_df)
    pssm_mat <- pssm_mat + prior
    pssm_mat <- pssm_mat / rowSums(pssm_mat)
    mat_to_pssm(pssm_mat)
}

#' Reverse complement a PSSM
#'
#' @param pssm A PSSM. Data frame with columns 'A', 'C', 'G', 'T' and 'pos' or a matrix with columns 'A', 'C', 'G', 'T'
#' @return A PSSM with the same format, but reverse complemented.
#'
#' @examples
#' # Create simulated PSSM data frame
#' pssm <- data.frame(
#'     pos = 1:4,
#'     A = c(0.1, 0.2, 0.3, 0.1),
#'     C = c(0.1, 0.3, 0.2, 0.1),
#'     G = c(0.1, 0.3, 0.3, 0.7),
#'     T = c(0.7, 0.2, 0.2, 0.1)
#' )
#'
#' # Reverse complement the PSSM
#' rc_pssm <- pssm_rc(pssm)
#'
#' @export
pssm_rc <- function(pssm) {
    pssm_orig <- pssm
    if (is.matrix(pssm_orig)) {
        pssm <- mat_to_pssm(pssm)
    }
    original_pos <- pssm$pos
    pssm <- pssm %>%
        mutate(tmp_A = A, tmp_C = C, tmp_G = G, tmp_T = T) %>%
        mutate(A = tmp_T, T = tmp_A, C = tmp_G, G = tmp_C) %>%
        select(-starts_with("tmp_")) %>%
        arrange(desc(pos)) %>%
        mutate(pos = original_pos)
    if (is.matrix(pssm_orig)) {
        pssm <- pssm_to_mat(pssm)
    }
    return(pssm)
}

#' Reverse Complement DNA Sequences
#'
#' This function takes a character vector of DNA sequences and returns their reverse complements.
#' It uses an efficient C++ implementation via Rcpp for improved performance.
#'
#' @param dna A character vector of DNA sequences. Can be a single sequence or multiple sequences.
#'            The sequences can be in upper or lower case.
#'
#' @return A character vector of the same length as the input, where each element
#'         is the reverse complement of the corresponding input sequence.
#'
#' @details The function performs the following operations on each sequence:
#'          1. Converts the sequence to uppercase.
#'          2. Reverses the sequence.
#'          3. Complements each base (A<->T, C<->G).
#'          Non-standard characters (not A, T, C, or G) are preserved in their reversed positions.
#'
#' @examples
#' rc("ATCG") # Returns "CGAT"
#' rc(c("ATCG", "GGCC", "TATA")) # Returns c("CGAT", "GGCC", "TATA")
#'
#' @export
rc <- function(dna) {
    if (!is.character(dna)) {
        cli_abort("The input should be a character vector")
    }
    rc_cpp(dna)
}

#' Trim PSSM
#'
#' This function trims a Position-Specific Scoring Matrix (PSSM) by removing positions with low information content at the beginning and end of the motif.
#'
#' @param pssm A data frame representing the PSSM, with columns for position (pos) and bits per position (bits).
#' @param bits_thresh The threshold value for bits per position. Positions with bits above this threshold will be kept, while positions with bits below this threshold at the beginning and the end of the motif will be removed. The default value is 0.1.
#'
#' @return A trimmed PSSM data frame, with positions filtered based on the bits threshold.
#'
#' @export
trim_pssm <- function(pssm, bits_thresh = 0.1) {
    pssm <- as.data.frame(pssm) %>%
        mutate(pos = 1:n() - 1)
    bits <- bits_per_pos(pssm)
    above_threshold <- bits > bits_thresh
    first_above <- which(above_threshold)[1] - 1
    last_above <- which(above_threshold)[length(which(above_threshold))] - 1
    pssm <- pssm %>%
        filter(pos >= first_above, pos <= last_above) %>%
        mutate(pos = 1:n() - 1)
    pssm
}

#' @export
#' @rdname trim_pssm
pssm_trim <- trim_pssm

#' Convert a dataset to a list of PSSM matrices
#'
#' @param dataset A data frame with columns 'motif', 'A', 'C', 'G', 'T'
#' @return A named list of matrices
#' @keywords internal
#' @noRd
dataset_to_pssm_list <- function(dataset) {
    if (!is.data.frame(dataset)) {
        cli::cli_abort("The {.field dataset} argument should be a data frame")
    }
    if (!all(c("A", "C", "G", "T", "motif") %in% colnames(dataset))) {
        cli::cli_abort("The {.field dataset} data frame should have columns {.val A}, {.val C}, {.val G}, {.val T}, and {.val motif}")
    }

    # Get unique motifs
    motifs <- unique(dataset$motif)

    # Convert dataset to list of matrices
    pssm_list <- lapply(motifs, function(m) {
        mat <- dataset[dataset$motif == m, c("A", "C", "G", "T")]
        as.matrix(mat)
    })

    names(pssm_list) <- motifs
    return(pssm_list)
}

#' Concatenate two PSSM matrices
#'
#' @param pssm1,pssm2 PSSM matrices
#' @param gap The gap between the two PSSM matrices, can be 0
#' @param trim Whether to trim the PSSM matrices
#' @param bits_thresh The threshold value for bits per position. Positions with bits above this threshold will be kept, while positions with bits below this threshold at the beginning and the end of the motif will be removed. The default value is 0.1.
#' @param orientation The orientation of the PSSMs. One of "ff" (forward-forward), "fr" (forward-reverse), "rf" (reverse-forward), "rr" (reverse-reverse). Default is "ff".
#'
#' @examples
#' # Basic concatenation of two motifs
#' pssm1 <- as.matrix(MOTIF_DB["HOMER.GATA3_2"])
#' pssm2 <- as.matrix(MOTIF_DB["JASPAR.CDX1"])
#' concat_motif <- pssm_concat(pssm1, pssm2)
#' plot_pssm_logo(concat_motif, title = "Concatenated GATA3 + CDX1")
#'
#' # Concatenation with a gap
#' concat_with_gap <- pssm_concat(pssm1, pssm2, gap = 3)
#' plot_pssm_logo(concat_with_gap, title = "GATA3 + CDX1 with 3bp gap")
#'
#' # Different orientations
#' concat_fr <- pssm_concat(pssm1, pssm2, orientation = "fr") # forward-reverse
#' plot_pssm_logo(concat_fr, title = "GATA3 forward + CDX1 reverse")
#' concat_rf <- pssm_concat(pssm1, pssm2, orientation = "rf") # reverse-forward
#' plot_pssm_logo(concat_rf, title = "GATA3 reverse + CDX1 forward")
#' concat_rr <- pssm_concat(pssm1, pssm2, orientation = "rr") # reverse-reverse
#' plot_pssm_logo(concat_rr, title = "Both motifs reversed")
#'
#' # Without trimming
#' concat_no_trim <- pssm_concat(pssm1, pssm2, trim = FALSE)
#' plot_pssm_logo(concat_no_trim, title = "Concatenated without trimming")
#'
#' @export
pssm_concat <- function(pssm1, pssm2, gap = 0, trim = TRUE, bits_thresh = 0.1, orientation = "ff") {
    if (is.matrix(pssm1)) {
        pssm1 <- mat_to_pssm(pssm1)
    }
    if (is.matrix(pssm2)) {
        pssm2 <- mat_to_pssm(pssm2)
    }

    # Apply reverse complement based on orientation
    valid_orientations <- c("ff", "fr", "rf", "rr")
    if (!orientation %in% valid_orientations) {
        cli_abort("Invalid orientation. Must be one of: {.val {valid_orientations}}")
    }

    # Apply reverse complement as needed
    if (orientation == "rf" || orientation == "rr") {
        pssm1 <- pssm_rc(pssm1)
    }

    if (orientation == "fr" || orientation == "rr") {
        pssm2 <- pssm_rc(pssm2)
    }

    if (trim) {
        pssm1 <- trim_pssm(pssm1, bits_thresh)
        pssm2 <- trim_pssm(pssm2, bits_thresh)
    }

    if (gap > 0) {
        gap_df <- data.frame(pos = 1:gap, A = 0.25, C = 0.25, G = 0.25, T = 0.25)
        pssm_concat <- dplyr::bind_rows(pssm1, gap_df, pssm2)
    } else {
        pssm_concat <- dplyr::bind_rows(pssm1, pssm2)
    }

    pssm_concat <- pssm_concat %>%
        mutate(pos = 1:n() - 1) %>%
        select(pos, A, C, G, T)

    return(pssm_concat)
}

#' @export
#' @rdname pssm_concat
concat_pssm <- pssm_concat

#' Screen sequences for positions with PWM scores meeting a threshold condition
#'
#' @description This function screens sequences to find positions where the local PWM score
#' meets a specified threshold condition using various operators (>, <, >=, <=, ==).
#'
#' @param operator A character string specifying the comparison operator. One of: ">", "<", ">=", "<=", "=="
#' @param threshold A numeric value specifying the threshold for comparison
#'
#' @return A list with one element per sequence, where each element is a numeric vector
#' containing the 1-indexed positions that meet the threshold condition.
#'
#' @examples
#' \dontrun{
#' # Find positions where HNF1A motif scores above -19
#' hnf1a_positions <- screen_local_pwm(
#'     cluster_sequences_example,
#'     as.matrix(MOTIF_DB["JASPAR.HNF1A"]),
#'     operator = ">",
#'     threshold = -37
#' )
#'
#' which(purrr::map_dbl(hnf1a_positions, length) > 0) # which sequences have positions?
#' hnf1a_positions[[4]] # positions for the 4th sequence
#'
#' hnf1a_energy <- compute_local_pwm(cluster_sequences_example, as.matrix(MOTIF_DB["JASPAR.HNF1A"]))
#' hnf1a_energy[4, hnf1a_positions[[4]]]
#' }
#'
#' @inheritParams compute_local_pwm
#' @export
screen_local_pwm <- function(sequences, pssm, operator, threshold, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01) {
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

    # Validate operator
    valid_operators <- c(">", "<", ">=", "<=", "==")
    if (!operator %in% valid_operators) {
        cli_abort("The {.field operator} must be one of: {.val {valid_operators}}")
    }

    pssm_mat <- as.matrix(pssm[, c("A", "C", "G", "T")])

    if (nrow(pssm_mat) == 0) {
        cli_abort("The {.field pssm} matrix should have at least one row")
    }

    if (prior < 0 || prior > 1) {
        cli_abort("The {.field prior} should be between 0 and 1")
    }

    if (prior > 0) {
        pssm_mat <- pssm_mat + prior
    }

    positions <- screen_local_pwm_cpp(
        sequences = toupper(sequences),
        pssm_mat = pssm_mat,
        is_bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_factor = spat$spat_factor,
        bin_size = binsize,
        operator_str = operator,
        threshold = threshold
    )

    return(positions)
}
