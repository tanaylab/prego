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

    if (is.null(spat_max)) {
        spat_max <- nchar(sequences[[1]])
    }

    if (is.null(spat_min)) {
        spat_min <- 1
    }

    if (!(spat_min == 1 && spat_max == nchar(sequences[[1]]))) {
        sequences <- stringr::str_sub(sequences, start = spat_min, end = spat_max)
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
#' @return a matrix with \code{length(sequences)} rows and \code{ncol(pssm)} columns with the local PWM for each sequence in each position.
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#'
#' pwm <- compute_local_pwm(cluster_sequences_example, res$pssm, res$spat)
#' head(pwm)
#' }
#'
#' @inheritParams compute_pwm
#' @export
compute_local_pwm <- function(sequences, pssm, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01) {
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

    seq_l <- stringr::str_length(sequences)
    if (length(unique(seq_l)) > 1) {
        cli_abort("All sequences should have the same length")
    }

    pwm <- compute_local_pwm_cpp(
        sequences = toupper(sequences),
        pssm_mat = pssm_mat,
        is_bidirect = bidirect,
        spat_min = spat_min,
        spat_max = spat_max,
        spat_factor = spat$spat_factor,
        bin_size = binsize
    )

    return(pwm)
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

pssm_mat_to_df <- function(pss_mat) {
    pssm_df <- as.data.frame(pss_mat)
    pssm_df$pos <- rownames(pssm_df)
    pssm_df <- pssm_df %>%
        select(pos, A, C, G, T)
    return(pssm_df)
}

pssm_add_prior <- function(pssm_df, prior) {
    pssm_mat <- pssm_to_mat(pssm_df)
    pssm_mat <- pssm_mat + prior
    pssm_mat <- pssm_mat / rowSums(pssm_mat)
    pssm_mat_to_df(pssm_mat)
}

#' Compute KL divergence between two PSSMs
#'
#' @param pssm1 first PSSM matrix or data frame
#' @param pssm2 second PSSM matrix or data frame
#'
#' @return KL divergence between the two PSSMs
#'
#' @examples
#' \dontrun{
#' res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' pssm_diff(res1$pssm, JASPAR_motifs[JASPAR_motifs$motif == "HNF1A", ])
#' }
#'
#' @export
pssm_diff <- function(pssm1, pssm2) {
    pssm1 <- pssm_to_mat(pssm1)
    pssm2 <- pssm_to_mat(pssm2)
    n_pos1 <- nrow(pssm1)
    n_pos2 <- nrow(pssm2)

    if (n_pos1 == 0 || n_pos2 == 0) {
        cli::cli_abort("PSSM matrices cannot be empty")
    }

    window_size <- min(n_pos1, n_pos2)
    max_pos <- max(n_pos1, n_pos2)


    if (n_pos1 > n_pos2) {
        pssm_s <- pssm2
        pssm_l <- pssm1
    } else {
        pssm_s <- pssm1
        pssm_l <- pssm2
    }

    epsilon <- 1e-10 # to avoid division by 0 or log(0)
    pssm_s <- pssm_s + epsilon
    pssm_s <- pssm_s / rowSums(pssm_s) # renormalize
    pssm_l <- pssm_l + epsilon
    pssm_l <- pssm_l / rowSums(pssm_l)

    kl <- function(a, b) {
        0.5 / ncol(a) * sum(colSums(a * log(a / b) +
            b * log(b / a)))
    }

    kl_scores <- purrr::map_dbl(1:(max_pos - window_size + 1), ~ {
        kl(t(pssm_l[.x:(.x + window_size - 1), ]), t(pssm_s))
    })

    return(min(kl_scores))
}

#' Compute the correlation between two given PSSMs
#'
#' @description The correlation is computed by shifting the shorter PSSM along the longer one
#' and computing the correlation at each position. The maximum correlation is returned.
#'
#' @param pssm1 first PSSM matrix or data frame
#' @param pssm2 second PSSM matrix or data frame
#' @param method method to use for computing the correlation. See \code{\link[stats]{cor}} for details.
#' @param prior a prior probability for each nucleotide.
#'
#' @return Correlation between the two PSSMs
#'
#' @examples
#' \dontrun{
#' res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' pssm_cor(res1$pssm, JASPAR_motifs[JASPAR_motifs$motif == "HNF1A", ])
#' }
#'
#' @export
pssm_cor <- function(pssm1, pssm2, method = "spearman", prior = 0.01) {
    pssm1 <- pssm_to_mat(pssm1)
    pssm2 <- pssm_to_mat(pssm2)
    n_pos1 <- nrow(pssm1)
    n_pos2 <- nrow(pssm2)

    if (n_pos1 == 0 || n_pos2 == 0) {
        cli::cli_abort("PSSM matrices cannot be empty")
    }

    window_size <- min(n_pos1, n_pos2)
    max_pos <- max(n_pos1, n_pos2)

    if (n_pos1 > n_pos2) {
        pssm_s <- pssm2
        pssm_l <- pssm1
    } else {
        pssm_s <- pssm1
        pssm_l <- pssm2
    }

    pssm_s <- pssm_s + prior
    pssm_s <- pssm_s / rowSums(pssm_s) # renormalize
    pssm_l <- pssm_l + prior
    pssm_l <- pssm_l / rowSums(pssm_l)

    scores <- purrr::map_dbl(1:(max_pos - window_size + 1), ~ {
        cor(as.vector(t(pssm_l[.x:(.x + window_size - 1), ])), as.vector(t(pssm_s)), method = method)
    })

    return(max(scores))
}

#' Compute a correlation matrix for a pssm dataset
#'
#' @param dataset a pssm dataset. A data frame with the columns 'motif', 'pos', 'A", 'C', 'G', 'T'
#' @param parallel whether to use parallel computing
#'
#' @return A correlation matrix for the PSSMs in the dataset
#'
#' @examples
#' \dontrun{
#' cm <- pssm_dataset_cor(JASPAM_motifs)
#' head(cm)
#' }
#'
#' @inheritParams pssm_cor
#' @export
pssm_dataset_cor <- function(dataset, method = "spearman", prior = 0.01, parallel = getOption("prego.parallel", TRUE)) {
    motifs <- unique(dataset$motif)
    motif_combs <- t(utils::combn(motifs, 2)) %>% as.data.frame()
    colnames(motif_combs) <- c("motif1", "motif2")

    pssm_cors <- safe_adply(motif_combs, 1, function(x) {
        tibble(cor = pssm_cor(
            dataset %>% filter(motif == x$motif1),
            dataset %>% filter(motif == x$motif2)
        ))
    },
    .parallel = parallel
    )

    # transform to a matrix while making sure all the
    pssm_mat <- pssm_cors %>%
        tidyr::complete(motif1 = motifs, motif2 = motifs, fill = list(cor = 0)) %>%
        tidyr::spread(motif2, cor) %>%
        tibble::column_to_rownames("motif1") %>%
        as.matrix()

    pssm_mat <- pssm_mat[motifs, motifs]


    # fill the lower triangle
    pssm_mat[lower.tri(pssm_mat)] <- t(pssm_mat)[lower.tri(pssm_mat)]

    # fill the diagonal
    diag(pssm_mat) <- 1

    return(pssm_mat)
}


#' Match PSSM to a directory of motifs
#'
#' @description Match a PSSM to a directory of motifs. The PSSM is matched to each motif in the directory by computing the correlation between the two PSSMs.
#'
#' @param pssm PSSM matrix or data frame
#' @param motifs a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name
#' @param best return the best match only
#' @param parallel use parallel processing. Set the number of cores using \code{set_parallel}.
#'
#' @return if \code{best} is \code{TRUE}, a string with the best match. Otherwise, a data frame with a row per motif and a column named 'cor' with its correlation to \code{pssm}. The data frame is sorted by descreasing correlation.
#'
#' @examples
#' \dontrun{
#' res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' head(pssm_match(res1$pssm, JASPAR_motifs))
#' pssm_match(res1$pssm, JASPAR_motifs, best = TRUE)
#' }
#'
#' @inheritParams pssm_cor
#'
#' @export
pssm_match <- function(pssm, motifs, best = FALSE, method = "spearman", parallel = getOption("prego.parallel", TRUE)) {
    if (!is.data.frame(motifs)) {
        cli_abort("The {.field motifs} argument should be a data frame")
    }
    if (!all(c("A", "C", "G", "T") %in% colnames(motifs))) {
        cli_abort("The {.field motifs} data frame should have columns {.val A}, {.val C}, {.val G}, {.val T}")
    }
    if (!("motif" %in% colnames(motifs))) {
        cli_abort("The {.field motifs} data frame should have a column {.val motif}")
    }

    res <- safe_ddply(motifs, "motif", function(x) {
        tibble(cor = pssm_cor(pssm, x, method = method))
    }, .parallel = parallel)

    res <- res %>% arrange(desc(cor))

    if (best) {
        return(res$motif[1])
    }

    return(res)
}

#' Plot LOGO of the pssm result from the regression
#'
#' @param pssm PSSM matrix or data frame
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param pos_bits_thresh Positions with bits above this threshold would be highlighted in red. If \code{NULL}, no positions would be highlighted.
#' @param revcomp whether to plot the reverse complement of the PSSM
#' @param method Height method, can be one of "bits" or "probability" (default:"bits")
#'
#' @return a ggplot object
#'
#' @examples
#' pssm <- data.frame(
#'     pos = seq(0, 9, by = 1),
#'     A = c(
#'         0.16252439817826936, 0.4519127838188067, 0, 1, 0, 0.9789171974522293,
#'         0.9743866100297978, 0.013113942843003835, 0.3734676916683981,
#'         0.32658771473191045
#'     ),
#'     C = c(
#'         0.43038386467143785, 0.13116231900388756, 0, 0, 0, 0, 0, 0.46975132995175056,
#'         0.1669956368169541, 0.29795679333680375
#'     ),
#'     G = c(
#'         0.22999349381912818, 0.002929742520705392, 1, 0, 0, 0, 0.012679896024852597,
#'         0.4808858097241123, 0.4248389777685435, 0.20458094742321709
#'     ),
#'     T = c(
#'         0.1770982433311646, 0.41399515465660036, 0, 0, 1, 0.0210828025477707,
#'         0.012933493945349648, 0.03624891748113324, 0.0346976937461043,
#'         0.17087454450806872
#'     )
#' )
#' plot_pssm_logo(pssm)
#' \dontrun{
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_pssm_logo(res$pssm)
#' }
#'
#' @export
plot_pssm_logo <- function(pssm, title = "Sequence model", subtitle = ggplot2::waiver(), pos_bits_thresh = NULL, revcomp = FALSE, method = "bits") {
    if (revcomp) {
        pssm <- pssm_rc(pssm)
    }
    pfm <- t(pssm_to_mat(pssm))
    p <- ggseqlogo::ggseqlogo(pfm, method = method) +
        ggtitle(title, subtitle = subtitle)
    if (!is.null(pos_bits_thresh)) {
        bits <- bits_per_pos(t(pfm))
        pos_mask <- bits > pos_bits_thresh
        rect_data <- tibble(
            x = which(pos_mask) - 0.5,
            xend = x + 1,
            y = 0,
            yend = max(bits)
        )
        p <- p + geom_rect(data = rect_data, aes(xmin = x, xmax = xend, ymin = y, ymax = yend), fill = "red", alpha = 0.1)
    }
    return(p)
}

#' Plot LOGO of pssm from dataset (e.g. "HOMER" or "JASPAR")
#'
#' @param motif the motif name (e.g. "GATA4")
#' @param dataset a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs}, \code{JASPAR_motifs} or all_motif_datasets()
#'
#' @return a ggplot object
#'
#' @examples
#'
#' plot_pssm_logo_dataset("JASPAR.Brachyury")
#'
#' plot_pssm_logo_dataset("GATA5", JASPAR_motifs)
#'
#' @inheritParams plot_pssm_logo
#' @export
plot_pssm_logo_dataset <- function(motif, dataset = all_motif_datasets(), title = motif, subtitle = ggplot2::waiver(), pos_bits_thresh = NULL, revcomp = FALSE, method = "bits") {
    motif_dataset <- dataset %>%
        filter(motif == !!motif)
    if (nrow(motif_dataset) == 0) {
        cli_abort("The motif {.val {motif}} was not found in the dataset")
    }
    plot_pssm_logo(motif_dataset, title = title, subtitle = subtitle, pos_bits_thresh = pos_bits_thresh, revcomp = revcomp, method = method)
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
    pssm <- pssm %>%
        mutate(tmp_A = A, tmp_C = C, tmp_G = G, tmp_T = T) %>%
        mutate(A = tmp_T, T = tmp_A, C = tmp_G, G = tmp_C) %>%
        select(-starts_with("tmp_")) %>%
        arrange(desc(pos)) %>%
        mutate(pos = n() - pos + 1)
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
