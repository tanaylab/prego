#' Extract pwm of intervals from a motif database
#'
#' @param intervals misha intervals set
#'
#' @return The intervals set with additional columns per motif, containing the pwm of each interval for each motif
#'
#' @inheritParams extract_pwm
#'
#' @examples
#' \dontrun{
#' library(misha)
#' gdb.init_examples()
#' pwms <- gextract_pwm("annotations")
#' pwms[, 1:20]
#' }
#'
#' @export
gextract_pwm <- function(intervals, motifs = NULL, dataset = all_motif_datasets(), spat = NULL, spat_min = 1, spat_max = NULL, bidirect = TRUE, prior = 0.01, func = "logSumExp", parallel = getOption("prego.parallel", TRUE)) {
    if (!requireNamespace("misha", quietly = TRUE)) {
        cli_abort("The {.field misha} package is required for this function. Please install it with {.code remotes::install_packages('tanaylab/misha')}.")
    }

    if (is.character(intervals)) {
        intervals <- misha::gintervals.load(intervals)
    }

    sequences <- misha::gseq.extract(intervals)

    res <- extract_pwm(sequences, motifs = motifs, dataset = dataset, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior, func = func, parallel = parallel)

    return(cbind(intervals, as.data.frame(res)))
}

#' Extract quantiles of pwm of intervals from a motif database
#'
#' @description Extract for each interval its quantile in the genome for each motif given its length. Note
#' that the quantiles are computed for each motif separately, and therefore this might be slow for intervals with un-normalized lengths.
#'
#' @param intervals misha intervals set
#' @param percision the percision of the quantiles. Default is 0.01, which means that the quantiles will be computed for every 1% of the pwm.
#'
#' @return a data frame with the quantiles of the pwm for each interval and motif. The quantiles columns would be of the form \{motif\}.q
#'
#' @examples
#' \dontrun{
#' library(misha)
#' gdb.init_examples()
#' gextract_pwm.quantile("annotations", motifs = c("JASPAR.CDX1", "JASPAR.CDX2"), dist_from_edge = 100)
#' }
#'
#' @inheritParams gextract_pwm
#' @inheritParams gpwm_quantiles
#'
#' @export
gextract_pwm.quantile <- function(intervals, motifs = NULL, dataset = all_motif_datasets(), percision = 0.01, spat = NULL, spat_min = 1, spat_max = NULL, bidirect = TRUE, prior = 0.01, func = "logSumExp", n_sequences = 1e4, dist_from_edge = 3e6, chromosomes = NULL, parallel = getOption("prego.parallel", TRUE)) {
    withr::local_options(gmax.data.size = 1e9)
    pwms <- gextract_pwm(intervals, motifs = motifs, dataset = dataset, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior, parallel = parallel)

    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    if (is.character(intervals)) {
        intervals <- misha::gintervals.load(intervals)
    }

    pwms <- pwms %>%
        mutate(l = end - start)
    breaks <- seq(0, 1, by = percision)

    if (length(unique(pwms$l)) > 1) {
        cli_alert_warning("The length of the intervals is not the same. The computation of the quantiles will be done for each length separately which might take a while.")
        cli_ul("Please consider using {.code misha.ext::gintervals.normalize} to normalize the intervals to the same length before running this function.")
    }

    n_motifs <- length(unique(dataset$motif))
    n_lengths <- length(unique(pwms$l))

    if (parallel) {
        if (n_motifs > n_lengths) {
            motif_p <- TRUE
            len_p <- FALSE
        } else {
            motif_p <- FALSE
            len_p <- TRUE
        }
    } else {
        motif_p <- FALSE
        len_p <- FALSE
    }

    # go over every length and compute quantiles
    pwm_q <- plyr::ddply(pwms, "l", function(pwms_l) {
        pwms_l <- tibble::as_tibble(pwms_l)

        res <- plyr::daply(dataset, "motif", function(x) {
            quantiles <- gpwm_quantiles(size = pwms_l$l[1], quantiles = breaks, pssm = x, n_sequences = n_sequences, dist_from_edge = dist_from_edge, chromosomes = chromosomes, func = func)
            val2quant <- approxfun(x = quantiles, y = breaks, rule = 2)
            val2quant(pwms_l[[x$motif[1]]])
        }, .parallel = motif_p)

        if (is.null(dim(res))) {
            res <- as.matrix(res)
            colnames(res) <- dataset$motif[1]
        } else {
            res <- as.matrix(as.data.frame(t(res)))
        }

        res
    }, .parallel = len_p)

    pwm_q <- pwm_q %>%
        select(-l)

    if (n_motifs == 1) {
        colnames(pwm_q) <- dataset$motif[1]
    }

    colnames(pwm_q) <- paste0(colnames(pwm_q), ".q")

    return(cbind(intervals, as.data.frame(pwm_q)))
}



#' Compute quantile of pwm for a given interval size
#'
#' @description Computes the quantile of the pwm for a given interval size by sampling random intervals from the genome. The number of sequences to sample can be specified with \code{n_sequences}.
#'
#' @param size size of the intervals to sample
#' @param pssm PSSM matrix or data frame
#' @param quantiles quantiles to compute. See \code{quantile} for more details.
#' @param n_sequences number of sequences to sample in order to compute the quantiles. The default is 1e4.
#'
#' @return a named vector with the quantiles of the pwm for the given interval size.
#'
#' @inheritParams compute_pwm
#' @inheritParams misha.ext::grandom_genome
#'
#' @examples
#' \dontrun{
#' library(misha)
#' library(dplyr)
#' gdb.init_examples()
#' pssm <- JASPAR_motifs %>%
#'     filter(motif == "JASPAR.CDX1") %>%
#'     select(-motif)
#' gpwm_quantiles(1000, seq(0, 1, 0.1), pssm, dist_from_edge = 100)
#' }
#'
#' @export
gpwm_quantiles <- function(size, quantiles, pssm, spat = NULL, spat_min = 1, spat_max = NULL, bidirect = TRUE, prior = 0.01, n_sequences = 1e4, dist_from_edge = 3e6, chromosomes = NULL, func = "logSumExp") {
    if (!requireNamespace("misha.ext", quietly = TRUE)) {
        cli_abort("The {.field misha.ext} package is required for this function. Please install it with {.code remotes::install_packages('tanaylab/misha.ext')}.")
    }

    r_intervals <- misha.ext::grandom_genome(size, n_sequences, dist_from_edge, chromosomes)
    sequences <- misha::gseq.extract(r_intervals)
    pwm <- compute_pwm(sequences, pssm, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior, func = func)

    return(quantile(pwm, probs = quantiles, na.rm = TRUE))
}
