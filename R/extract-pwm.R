#' @rdname extract_pwm
#' @export
extract_pwm_old <- function(sequences, motifs = NULL, dataset = all_motif_datasets(), spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01, func = "logSumExp", parallel = getOption("prego.parallel", TRUE)) {
    if (inherits(dataset, "MotifDB")) {
        dataset <- as.data.frame(dataset)
    }

    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    sequences <- toupper(sequences)

    res <- safe_daply(dataset, "motif", function(x) {
        if ("motif" %in% colnames(spat)) {
            spat <- spat %>% filter(motif == x$motif[1])
        }
        compute_pwm(sequences, x, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior, func = func)
    }, .parallel = parallel)

    if (is.null(dim(res))) {
        res <- as.matrix(res)
        colnames(res) <- dataset$motif[1]
    } else {
        res <- as.matrix(as.data.frame(t(res)))
    }

    if (!is.null(names(sequences))) {
        rownames(res) <- names(sequences)
    }

    return(res)
}

#' Extract pwm of sequences from a motif database
#'
#' @description Extracts the pwm of a motif from a motif database. \code{extract_pwm_old} is a deprecated version of this function, which is slower, and returns slightly different results due to float percision instead of double. If the sequences are not of the same length, the old version will be used.
#'
#' @param motifs names of specific motifs to extract from the dataset
#' @param dataset a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs} or \code{JASPAR_motifs}, or \code{all_motif_datasets()}, or a MotifDB object.
#' @param parallel logical, whether to use parallel processing
#'
#' @return a matrix size of # of sequences x # of motifs with the pwm of each sequence for each motif
#'
#' @examples
#' \dontrun{
#' pwms <- extract_pwm(
#'     cluster_sequences_example,
#'     motifs = c("JASPAR.CDX1", "HOMER.Hnf1", "HOMER.GATA3_2")
#' )
#' head(pwms)
#'
#' # all motifs
#' all_pwms <- extract_pwm(cluster_sequences_example, prior = 0.01)
#' dim(all_pwms)
#' all_pwms[1:5, 1:5]
#'
#' # for a specific dataset
#' pwms_jaspar <- extract_pwm(cluster_sequences_example, dataset = JASPAR_motifs, prior = 0.01)
#' head(pwms_jaspar)
#'
#' # for specific motifs
#' pwms_jaspar <- extract_pwm(
#'     cluster_sequences_example,
#'     motifs = c("JASPAR.CDX1", "JASPAR.CDX2"),
#'     prior = 0.01
#' )
#' }
#'
#' @inheritParams compute_pwm
#' @export
extract_pwm <- function(sequences, motifs = NULL, dataset = MOTIF_DB, spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01, func = "logSumExp", parallel = getOption("prego.parallel", TRUE)) {
    # if not all sequences have the same length or the function is max, use the old version
    if (length(unique(nchar(sequences))) != 1 || func == "max") {
        return(extract_pwm_old(sequences, motifs = motifs, dataset = dataset, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, prior = prior, func = func, parallel = parallel))
    }

    if (inherits(dataset, "MotifDB")) {
        mdb <- dataset
        if (!is.null(motifs)) {
            mdb <- mdb[motifs]
        }
        if (mdb@prior != prior) {
            prior(mdb) <- prior
        }
    } else {
        if (!is.null(motifs)) {
            dataset <- dataset %>% filter(motif %in% motifs)
        }
        mdb <- create_motif_db(dataset, prior = prior)
    }

    if (!is.null(spat)) {
        mdb@spat_factors <- spat$spat_factor
        mdb@spat_bin_size <- unique(diff(spat$bin))[1]
    }

    if (!parallel) {
        RcppParallel::setThreadOptions(numThreads = 1)
        n_threads <- getOption("prego.parallel.nc", 1)
        withr::defer(RcppParallel::setThreadOptions(numThreads = n_threads))
    }

    if (is.null(spat_max) || is.na(spat_max)) {
        spat_max <- nchar(sequences[[1]])
    }

    if (is.null(spat_min) || is.na(spat_min)) {
        spat_min <- 1
    }

    if (!(spat_min == 1 && spat_max == nchar(sequences[[1]]))) {
        sequences <- stringr::str_sub(sequences, start = spat_min, end = spat_max)
    }

    res <- calc_seq_pwm(sequences, mdb, bidirect = bidirect)

    if (is.null(dim(res))) {
        res <- as.matrix(res)
        colnames(res) <- dataset$motif[1]
    }

    if (!is.null(names(sequences))) {
        rownames(res) <- names(sequences)
    }

    return(res)
}
