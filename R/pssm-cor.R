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
    pssm1 <- pssm_mat_to_df(pssm_to_mat(pssm1)) %>%
        mutate(motif = "pssm1") %>%
        select(motif, everything())
    pssm2 <- pssm_mat_to_df(pssm_to_mat(pssm2)) %>%
        mutate(motif = "pssm2") %>%
        select(motif, everything())

    if (nrow(pssm1) == 0 || nrow(pssm2) == 0) {
        cli::cli_abort("PSSM matrices cannot be empty")
    }

    return(pssm_dataset_cor(pssm1, pssm2, method = method, prior = prior)[1])
}

#' Compute correlation matrix for PSSM datasets using parallel processing
#'
#' @param dataset1 First PSSM dataset. A data frame with columns 'motif', 'A', 'C', 'G', 'T'
#' @param dataset2 Optional second PSSM dataset with the same structure as dataset1.
#'                If provided, computes correlations between motifs in dataset1 and dataset2.
#'                If NULL (default), computes correlations between all motifs in dataset1.
#' @param method Method to use for computing correlation. Either "spearman" (default) or "pearson"
#' @param prior A prior probability added to each nucleotide frequency (default: 0.01)
#'
#' @return If dataset2 is NULL, returns a symmetric square matrix of correlations between
#'         all motifs in dataset1. If dataset2 is provided, returns a matrix with
#'         rows corresponding to motifs in dataset1 and columns to motifs in dataset2.
#'
#' @examples
#' \dontrun{
#' # Correlations within a single dataset
#' cm <- pssm_dataset_cor(JASPAR_motifs)
#'
#' # Correlations between two datasets
#' cm2 <- pssm_dataset_cor(JASPAR_motifs, HOMER_motifs)
#' }
#'
#' @export
pssm_dataset_cor <- function(dataset1, dataset2 = NULL, method = c("spearman", "pearson"), prior = 0.01) {
    method <- match.arg(method)

    # Convert first dataset to list
    pssm_list1 <- dataset_to_pssm_list(dataset1)
    motifs1 <- names(pssm_list1)

    if (is.null(dataset2)) {
        # Compute correlations within dataset1
        cor_mat <- pssm_dataset_cor_parallel(pssm_list1, method = method, prior = prior)
        rownames(cor_mat) <- motifs1
        colnames(cor_mat) <- motifs1
    } else {
        # Convert second dataset to list
        pssm_list2 <- dataset_to_pssm_list(dataset2)
        motifs2 <- names(pssm_list2)

        # Compute correlations between datasets
        # We'll create a combined list but keep track of which motifs belong to which dataset
        combined_list <- c(pssm_list1, pssm_list2)
        n1 <- length(pssm_list1)
        n2 <- length(pssm_list2)

        cor_mat_full <- pssm_dataset_cor_parallel(combined_list, method = method, prior = prior)

        # Extract the relevant submatrix (dataset1 vs dataset2)
        cor_mat <- cor_mat_full[1:n1, (n1 + 1):(n1 + n2), drop = FALSE]
        rownames(cor_mat) <- motifs1
        colnames(cor_mat) <- motifs2
    }

    return(cor_mat)
}

#' Match PSSM to a directory of motifs using parallel processing
#'
#' @description Match a PSSM to a directory of motifs. The PSSM is matched to each motif
#' in the directory by computing the correlation between the two PSSMs.
#'
#' @param pssm PSSM matrix or data frame with columns 'A', 'C', 'G', 'T'
#' @param motifs A data frame with PSSMs ('A', 'C', 'G', 'T' columns), with an additional
#'        column 'motif' containing the motif name
#' @param best Whether to return only the best match (default: FALSE)
#' @param method Method to use for correlation calculation ("spearman" or "pearson")
#' @param prior A prior probability added to each nucleotide frequency
#'
#' @return If best is TRUE, returns a string with the best matching motif name.
#'         If best is FALSE, returns a data frame with columns 'motif' and 'cor',
#'         sorted by decreasing correlation.
#'
#' @examples
#' \dontrun{
#' # Find all matches
#' matches <- pssm_match(my_pssm, JASPAR_motifs)
#' head(matches)
#'
#' # Find best match only
#' best_match <- pssm_match(my_pssm, JASPAR_motifs, best = TRUE)
#' }
#'
#' @export
pssm_match <- function(pssm, motifs, best = FALSE, method = c("spearman", "pearson"), prior = 0.01) {
    method <- match.arg(method)

    # Convert pssm to proper format if needed
    if (is.data.frame(pssm)) {
        if (!all(c("A", "C", "G", "T") %in% colnames(pssm))) {
            cli::cli_abort("The {.field pssm} data frame should have columns {.val A}, {.val C}, {.val G}, {.val T}")
        }
        pssm <- as.matrix(pssm[, c("A", "C", "G", "T")])
    }

    # Create a temporary data frame with the query PSSM
    query_df <- pssm_mat_to_df(pssm) %>%
        mutate(motif = "query") %>%
        select(motif, everything())

    # Compute correlations between query and motif database
    cor_mat <- pssm_dataset_cor(query_df, motifs, method = method, prior = prior)

    # Extract results
    results <- data.frame(
        motif = colnames(cor_mat),
        cor = cor_mat[1, ],
        stringsAsFactors = FALSE
    )

    # Sort by correlation in descending order
    results <- results[order(results$cor, decreasing = TRUE), ]
    rownames(results) <- NULL

    if (best) {
        return(results$motif[1])
    }

    return(results)
}
