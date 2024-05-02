#' Get a data frame of all the motif datasets bundled with prego
#'
#' @description The data frame contain the PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name. Individual datasets are available within the package as \code{HOMER_motifs}, \code{JASPAR_motifs}, \code{JOLMA_motifs}, and \code{HOCOMOCO_motifs}.
#'
#' @return a data frame which concatenates motifs from "HOMER", "JASPAR" and "JOLMA".
#' Motif names are prefixed with the dataset name, e.g. "JASPAR.GATA4".
#'
#' @examples
#' all_motif_datasets()
#'
#' @references
#' \describe{
#' \item{HOMER: }{Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432}
#' \item{JASPAR: }{Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, Berhanu Lemma R, Turchi L, Blanc-Mathieu R, Lucas J, Boddie P, Khan A, Manosalva Pérez N, Fornes O, Leung TY, Aguirre A, Hammal F, Schmelter D, Baranasic D, Ballester B, Sandelin A, Lenhard B, Vandepoele K, Wasserman WW, Parcy F, and Mathelier A JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles Nucleic Acids Res. 2022 Jan 7;50(D1):D165-D173.; doi: 10.1093/nar/gkab1113}
#' \item{JOLMA: }{Jolma, A., Yin, Y., Nitta, K. et al. DNA-dependent formation of transcription factor pairs alters their binding specificity. Nature 534, S15–S16 (2016). \url{https://doi.org/10.1038/nature18912}}
#' \item{HOCOMOCO: }{Ivan V. Kulakovskiy; Ilya E. Vorontsov; Ivan S. Yevshin; Ruslan N. Sharipov; Alla D. Fedorova; Eugene I. Rumynskiy; Yulia A. Medvedeva; Arturo Magana-Mora; Vladimir B. Bajic; Dmitry A. Papatsenko; Fedor A. Kolpakov; Vsevolod J. Makeev: HOCOMOCO: towards a complete collection of transcription factor binding models for human and mouse via large-scale ChIP-Seq analysis. Nucl. Acids Res., Database issue, gkx1106 (11 November 2017). \url{https://doi.org/10.1093/nar/gkx1106}}
#' }
#'
#'
#' @export
all_motif_datasets <- function() {
    dplyr::bind_rows(
        prego::HOMER_motifs %>% mutate(dataset = "HOMER"),
        prego::JASPAR_motifs %>% mutate(dataset = "JASPAR"),
        prego::JOLMA_motifs %>% mutate(dataset = "JOLMA"),
        prego::HOCOMOCO_motifs %>% mutate(dataset = "HOCOMOCO")
    ) %>%
        mutate(motif_orig = motif) %>%
        tidyr::unite("motif", dataset, motif, sep = ".", remove = FALSE)
}

#' PSSMs from the HOMER motif database
#'
#'
#' @format A data frame containing the PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name.
#'
#' @references Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
#'
#' @source \url{http://homer.ucsd.edu/homer/motif/}
"HOMER_motifs"

#' PSSMs from the JASPAR motif database
#'
#'
#' @format A data frame containing the PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name.
#'
#' @references Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, Berhanu Lemma R, Turchi L, Blanc-Mathieu R, Lucas J, Boddie P, Khan A, Manosalva Pérez N, Fornes O, Leung TY, Aguirre A, Hammal F, Schmelter D, Baranasic D, Ballester B, Sandelin A, Lenhard B, Vandepoele K, Wasserman WW, Parcy F, and Mathelier A JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles Nucleic Acids Res. 2022 Jan 7;50(D1):D165-D173.; doi: 10.1093/nar/gkab1113
#'
#' @source \url{https://jaspar.genereg.net/downloads/}
"JASPAR_motifs"

#' PSSMs from the Jolma et al. motif database
#'
#'
#' @format A data frame containing the PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name.
#'
#' @references Jolma, A., Yin, Y., Nitta, K. et al. DNA-dependent formation of transcription factor pairs alters their binding specificity. Nature 534, S15–S16 (2016). \url{https://doi.org/10.1038/nature18912}
#'
#' @source \url{https://doi.org/10.1038/nature18912}
"JOLMA_motifs"

#' PSSMs from the HOCOMOCO motif database
#'
#'
#' @format A data frame containing the PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name.
#'
#' @references Ivan V. Kulakovskiy; Ilya E. Vorontsov; Ivan S. Yevshin; Ruslan N. Sharipov; Alla D. Fedorova; Eugene I. Rumynskiy; Yulia A. Medvedeva; Arturo Magana-Mora; Vladimir B. Bajic; Dmitry A. Papatsenko; Fedor A. Kolpakov; Vsevolod J. Makeev: HOCOMOCO: towards a complete collection of transcription factor binding models for human and mouse via large-scale ChIP-Seq analysis. Nucl. Acids Res., Database issue, gkx1106 (11 November 2017). \url{https://doi.org/10.1093/nar/gkx1106}
#'
#' @source \url{https://hocomoco11.autosome.ru/downloads/}
"HOCOMOCO_motifs"

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
extract_pwm <- function(sequences, motifs = NULL, dataset = all_motif_datasets(), spat = NULL, spat_min = 0, spat_max = NULL, bidirect = TRUE, prior = 0.01, func = "logSumExp", parallel = getOption("prego.parallel", TRUE)) {
    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    sequences <- toupper(sequences)

    res <- plyr::daply(dataset, "motif", function(x) {
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



#' Extract pssm of sequences from a motif database
#'
#' @param motif name of the motif to extract from the dataset
#' @param dataset a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs} or \code{JASPAR_motifs}, or \code{all_motif_datasets()}.
#'
#' @return a data frame with the pssm of the motif
#'
#' @examples
#' get_motif_pssm("JASPAR.HNF1A")
#'
#' @export
get_motif_pssm <- function(motif, dataset = all_motif_datasets()) {
    dataset %>%
        filter(motif == !!motif) %>%
        select(pos, A, C, G, T)
}
