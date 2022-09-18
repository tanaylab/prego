#' Get a data frame of all the motif datasets bundled with prego
#'
#' @description The data frame contain the PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name. Individual datasets are available within the package as \code{HOMER_motifs}, \code{JASPAR_motifs}, and \code{JOLMA_motifs}.
#'
#' @return a data frame which concatenates motifs from "HOMER", "JASPAR" and "JOLMA".
#' Motif names are prefixed with the dataset name, e.g. "JASPAR.GATA4".
#'
#' @examples
#' all_motif_datasets()
#'
#' @references
#' \itemize{
#' \item{HOMER: } {Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432}
#' \item{JASPAR: }{Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, Berhanu Lemma R, Turchi L, Blanc-Mathieu R, Lucas J, Boddie P, Khan A, Manosalva Pérez N, Fornes O, Leung TY, Aguirre A, Hammal F, Schmelter D, Baranasic D, Ballester B, Sandelin A, Lenhard B, Vandepoele K, Wasserman WW, Parcy F, and Mathelier A JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles Nucleic Acids Res. 2022 Jan 7;50(D1):D165-D173.; doi: \link{10.1093/nar/gkab1113}}
#' \item{JOLMA: }{Jolma, A., Yin, Y., Nitta, K. et al. DNA-dependent formation of transcription factor pairs alters their binding specificity. Nature 534, S15–S16 (2016). \link{https://doi.org/10.1038/nature18912}}
#' }
#'
#'
#' @export
all_motif_datasets <- function() {
    dplyr::bind_rows(
        HOMER_motifs %>% mutate(dataset = "HOMER"),
        JASPAR_motifs %>% mutate(dataset = "JASPAR"),
        JOLMA_motifs %>% mutate(dataset = "JOLMA")
    ) %>%
        mutate(motif_orig = motif) %>%
        tidyr::unite("motif", dataset, motif, sep = ".", remove = FALSE)
}

#' Extract pwm of sequences from a motif database
#'
#' @description Extracts the pwm of a motif from a motif database. Note that for a large number of motifs, this function can be slow and consume
#' a lot of memory.
#'
#' @param motifs names of specific motifs to extract from the dataset
#' @param dataset a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs} or \code{JASPAR_motifs}, or \code{all_motif_datasets()}.
#' @param parallel logical, whether to use parallel processing
#'
#' @return a matrix size of # of sequences x # of motifs with the pwm of each sequence for each motif
#'
#' @examples
#' pwms <- extract_pwm(cluster_sequences_example, motifs = c("JASPAR.CDX1", "HOMER.Hnf1", "HOMER.GATA3_2"))
#' head(pwms)
#'
#' # all motifs
#' all_pwms <- extract_pwm(cluster_sequences_example)
#' dim(all_pwms)
#' all_pwms[1:5, 1:5]
#'
#' # for a specific dataset
#' pwms_jaspar <- extract_pwm(cluster_sequences_example, dataset = JASPAR_motifs)
#' head(pwms_jaspar)
#'
#' @inheritParams compute_pwm
#' @export
extract_pwm <- function(sequences, motifs = NULL, dataset = all_motif_datasets(), spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, parallel = getOption("prego.parallel", TRUE)) {
    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    sequences <- toupper(sequences)

    res <- plyr::daply(dataset, "motif", function(x) {
        compute_pwm(sequences, x, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect)
    }, .parallel = parallel)

    res <- as.matrix(t(res))

    if (!is.null(names(sequences))) {
        rownames(res) <- names(sequences)
    }

    return(res)
}

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
gextract_pwm <- function(intervals, motifs = NULL, dataset = all_motif_datasets(), spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, parallel = getOption("prego.parallel", TRUE)) {
    if (!("misha" %in% installed.packages())) {
        cli_abort("The {.field misha} package is required for this function. Please install it with {.code remotes::install_packages('tanaylab/misha')}.")
    }

    if (is.character(intervals)) {
        intervals <- misha::gintervals.load(intervals)
    }

    sequences <- misha::gseq.extract(intervals)

    res <- extract_pwm(sequences, motifs = motifs, dataset = dataset, spat = spat, spat_min = spat_min, spat_max = spat_max, bidirect = bidirect, parallel = parallel)

    return(cbind(intervals, as.data.frame(res)))
}
