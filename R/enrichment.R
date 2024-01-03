#' Calculate motif enrichment
#'
#' Calculates motif enrichment for groups of loci (usually clusters).
#'
#' @param pwm_q A matrix representing the quantile values of a position weight matrix (PWM). This should be the output of \code{\link{gextract_pwm.quantile}}.
#' @param groups A character vector or factor specifying the groups to which each sequence belongs. Usually a clustering result.
#' @param threshold The threshold value for considering a position in the PWM as significant. Default is 0.99.
#' @param type The type of enrichment calculation to perform. Possible values are "relative" (default) or "absolute".
#'
#' @return A matrix representing the motif enrichment values for each group. When \code{type = "relative"}, the values are the relative enrichment of the motif in the group compared the loci in all the other groups. When \code{type = "absolute"}, the values are the enrichment of the motif compared to random genome.
#'
#'
#' @examples
#' \dontrun{
#' library(misha)
#' gdb.init_examples()
#' annot <- misha.ext::gintervals.normalize(gintervals.load("annotations"), 300)
#' pwm_q <- gextract_pwm.quantile(annot, motifs = c("JASPAR.CDX1", "JASPAR.CDX2"), dist_from_edge = 100)
#' pwm_q <- as.matrix(pwm_q[, c("JASPAR.CDX1.q", "JASPAR.CDX2.q")])
#' groups <- c("Group1", "Group1", "Group2", "Group2", "Group2", "Group3", "Group3", "Group3")
#'
#' # The threshold of 0.1 is used for demonstration purposes only. In practice, a threshold of 0.99 is recommended.
#' motif_enrichment(pwm_q, groups, threshold = 0.1, type = "relative")
#' motif_enrichment(pwm_q, groups, threshold = 0.1, type = "absolute")
#' }
#'
#' @export
motif_enrichment <- function(pwm_q, groups, threshold = 0.99, type = "relative") {
    if (!is.matrix(pwm_q)) {
        cli_abort("{.field pwm_q} must be a matrix")
    }

    if (length(groups) != nrow(pwm_q)) {
        cli_abort("The number of rows in {.field pwm_q} must be the same as the length of {.field groups}")
    }

    if (!is.character(groups) && !is.factor(groups)) {
        cli_abort("{.field groups} must be a character vector or a factor")
    }

    pwm_mf <- pwm_q >= threshold
    group_motifs <- tgs_matrix_tapply(t(pwm_mf + 0), groups, sum)
    group_n <- table(groups)[rownames(group_motifs)]

    # n_fg_ok
    n_fg_ok <- group_motifs

    # n_fg is group_n replicated to match the dimensions of group_motifs
    n_fg <- matrix(rep(group_n, each = ncol(group_motifs)), nrow = nrow(group_motifs), ncol = ncol(group_motifs), byrow = TRUE)

    # Total occurrences of each motif across all groups
    total_motif_occurrences <- colSums(group_motifs)

    # n_bg_ok is the total occurrences minus the occurrences in the current group
    n_bg_ok <- matrix(rep(total_motif_occurrences, nrow(group_motifs)), nrow = nrow(group_motifs), byrow = TRUE) - n_fg_ok

    # Total loci in all groups
    total_loci <- sum(group_n)

    # n_bg is the total loci minus the loci in the current group
    n_bg <- matrix(rep(total_loci, nrow(group_motifs)), nrow = nrow(group_motifs), ncol = ncol(group_motifs), byrow = TRUE) - n_fg

    # Calculate enrichment
    if (type == "relative") {
        rel_enrichment <- (n_fg_ok / n_fg) / (n_bg_ok / n_bg)
        return(rel_enrichment)
    } else if (type == "absolute") {
        abs_enrcihment <- (n_fg_ok / n_fg) / (1 - threshold)
        return(abs_enrcihment)
    } else {
        cli_abort("{.field type} must be one of {.val relative} or {.val absolute}")
    }
}
