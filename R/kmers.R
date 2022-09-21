
#' Screen for kmers
#'
#' @param sequences A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through \code{toupper})
#' @param response A matrix of response variables - number of rows should equal the number of sequences
#' @param kmer_length The number of non-gap characters in motifs that will be screened
#' @param min_cor Only patterns for which the maximum correlation to one of the response variable is larger than min_cor will be reported
#' @param min_n Only patterns for which the average number of occurrences in the sequences is larger than min_n will be reported
#' @param is_train a boolean vector that determine which subset of sequences to use when screening
#' @param min_gap,max_gap the length of a gap to be considered in the pattern. Only one gap, of length min_gap:max_gap, is being used, and is located anywhere in the motif. Note that this greatly expand the search space (and increase multiple testing severely).
#' @param from_range Sequences will be considered only from position from_range (default 0)
#' @param to_range Sequences will be considered only up to position to_range (default NULL - using the length of the sequences)
#' @param return_mat Return a matrix of patterns and their correlation to the response variables instead of a data frame. (default: FALSE)
#' @param seed random seed
#' @param verbose show verbose messages
#'
#' @return A data frame with the following columns, together with a column for each response variable with the
#' correlation of the kmers to the response variable:
#' \itemize{
#' \item{kmer: }{the kmer pattern, where "*" indicates a wildcard},
#' \item{max_r2: }{the maximum R^2 to one of the response variables},
#' \item{avg_n: }{the average number of times the kmer appears in the sequences},
#' \item{avg_var: }{the variance of the number of times the kmer appears in the sequences},
#' }
#' if \code{return_mat} is TRUE, a matrix with correlations to the response variables (where
#' rows are the kmers) is returned instead of a data frame.
#' If no kmer is found, an empty data frame is returned.
#'
#' @examples
#' kmers <- screen_kmers(sequences_example, response_mat_example)
#' head(kmers)
#'
#' kmers <- screen_kmers(sequences_example, response_mat_example, return_mat = TRUE)
#' head(kmers)
#'
#' kmers <- screen_kmers(sequences_example, response_mat_example, max_gap = 3)
#' head(kmers)
#'
#' @export
screen_kmers <- function(sequences,
                         response,
                         kmer_length = 6,
                         min_cor = 0.08,
                         min_n = 50,
                         is_train = NULL,
                         min_gap = 0,
                         max_gap = 0,
                         from_range = 0,
                         to_range = NULL,
                         return_mat = FALSE,
                         seed = 60427,
                         verbose = FALSE) {
    if (is.null(is_train)) {
        is_train <- rep(TRUE, length(sequences))
    }

    if (length(sequences) != length(is_train)) {
        cli_abort("The number of sequences and the length of {.field is_train} vector do not match")
    }

    if (is.null(nrow(response))) {
        response <- matrix(response, ncol = 1)
    }

    if (length(sequences) != nrow(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (length(kmer_length) > 1) {
        cli_abort("The {.field kmer_length} argument should be a single number. If you are want to run multiple kmer lengths, please set the {.field multi_kmers} to TRUE.")
    }

    n_in_train <- sum(is_train)

    if (is.null(to_range)) {
        to_range <- nchar(sequences[1])
    }

    cli_alert_info("Number of response variables: {.val {ncol(response)}}")
    cli_alert_info("Screening kmers of length {.val {kmer_length}}, from position {.val {from_range}} to position {.val {to_range}}")
    if (max_gap < min_gap) {
        cli_abort("{.field min_gap} should be smaller than {.field max_gap}")
    }

    if (min_gap < 0) {
        cli_abort("{.field min_gap} should be larger than 0")
    }

    if (!all(min_gap:max_gap == 0)) {
        cli_alert_info("Gaps of length {.val {min_gap}}:{.val {max_gap}} are allowed")
    }


    cli_alert_info("minimal correlation: {.val {min_cor}}, minimal number of occurrences: {.val {min_n}}")

    res <- screen_kmers_cpp(
        toupper(sequences),
        response,
        is_train,
        kmer_length,
        from_range,
        to_range,
        min_cor,
        min_n,
        min_gap,
        max_gap,
        n_in_train,
        seed,
        verbose
    )

    cli_alert_success("Found {.val {nrow(res)}} kmers in {.val {length(sequences)}} sequences.")

    if (!is.data.frame(res)) {
        cli_warn("did not find any kmers")
        res <- as.data.frame(res)
    }

    res <- res %>%
        mutate(kmer = replace_gaps_with_wildcards(kmer))

    res <- res %>%
        arrange(desc(max_r2))

    if (return_mat) {
        res_mat <- as.matrix(res[, -(1:4)])
        rownames(res_mat) <- res$kmer
        return(res_mat)
    }

    return(res)
}

replace_gaps_with_wildcards <- function(seqs, gap_char = "*") {
    stringr::str_replace(seqs, "\\d+", function(x) {
        purrr::map_chr(x, function(s) {
            if (is.na(s)) {
                ""
            } else {
                paste(rep(gap_char, s), collapse = "")
            }
        })
    })
}
