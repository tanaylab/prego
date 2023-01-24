.onLoad <- function(libname, pkgname) {
    limit_cores <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(limit_cores) && limit_cores == "TRUE") {
        # use 2 cores in CRAN
        set_parallel(2L)
    } else {
        set_parallel()
    }

    # this is here in order to avoid R CMD check warnings
    globals <- c("motif", "pos", "A", "C", "G", "T", "HOMER_motifs", "JASPAR_motifs", "JOLMA_motifs", "dataset", "max_r2", "kmer", "kmer_clust", "idx", "fold", "resp", "score", "comb_score", "model", "label", "bin", "spat_factor", "i", "pred", "dist", "pos", ".", "cluster", "sequences", "value", "name", "cluster_ids", "clusters", "y")
    utils::suppressForeignCheck(globals)
    utils::globalVariables(globals)
}
