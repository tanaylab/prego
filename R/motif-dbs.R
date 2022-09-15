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
