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
