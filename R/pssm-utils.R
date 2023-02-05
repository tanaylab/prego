#' Compute PWMs for a set of sequences given a PSSM matrix
#'
#' @param sequences a vector of sequences
#' @param pssm a PSSM matrix or data frame. The columns of the matrix or data frame should be named with the nucleotides ('A', 'C', 'G' and 'T').
#' @param spat a data frame with the spatial model (as returned from the \code{$spat} slot from the regression). Should contain a column called 'bin' and a column called 'spat_factor'.
#' @param bidirect is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.
#' @param prior a prior probability for each nucleotide.
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
#' @inheritParams regress_pwm
#' @export
compute_pwm <- function(sequences, pssm, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0) {
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

    pwm <- compute_pwm_cpp(
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
consensus_from_pssm <- function(pssm, single_thresh = 0.6, double_thresh = 0.85) {
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

    pssm_df %>%
        arrange(as.numeric(pos)) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("pos") %>%
        select(A, C, G, T) %>%
        as.matrix()
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
pssm_cor <- function(pssm1, pssm2, method = "spearman") {
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

    scores <- purrr::map_dbl(1:(max_pos - window_size + 1), ~ {
        cor(as.vector(t(pssm_l[.x:(.x + window_size - 1), ])), as.vector(t(pssm_s)), method = method)
    })

    return(max(scores))
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

    res <- plyr::ddply(motifs, "motif", function(x) {
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
#' @param pssm the 'pssm' field from the regression result
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#'
#' @return a ggplot object
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(sequences_example, response_mat_example)
#' plot_pssm_logo(res$pssm)
#' }
#'
#' @export
plot_pssm_logo <- function(pssm, title = "Sequence model", subtitle = ggplot2::waiver()) {
    pfm <- t(pssm_to_mat(pssm))
    ggseqlogo::ggseqlogo(pfm) +
        ggtitle(title, subtitle = subtitle)
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
plot_pssm_logo_dataset <- function(motif, dataset = all_motif_datasets(), title = motif, subtitle = ggplot2::waiver()) {
    motif_dataset <- dataset %>%
        filter(motif == !!motif)
    if (nrow(motif_dataset) == 0) {
        cli_abort("The motif {.val {motif}} was not found in the dataset")
    }
    plot_pssm_logo(motif_dataset, title = title, subtitle = subtitle)
}
