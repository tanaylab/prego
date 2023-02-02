#' Screen for motifs in a database given a response variable
#'
#' @param sequences a vector with the sequences
#' @param response a vector of response variable for each sequence. If the response is a matrix, the average will be used.
#' @param metric metric to use in order to choose the best motif. One of 'ks' or 'r2'. If NULL - the default would be 'ks' for binary variables, and 'r2' for continuous variables.
#' @param only_best return only the best motif (the one with the highest score). If FALSE, all the motifs will be returned.
#' @param alternative alternative hypothesis for the KS test. One of 'two.sided', 'less' or 'greater'.
#'
#'
#' @inheritParams extract_pwm
#' @inheritDotParams compute_pwm
#'
#' @return a data frame with the following columns:
#' \itemize{
#' \item{motif: }{the motif name.}
#' \item{score: }{the score of the motif (depending on \code{metric}).}
#' }
#'
#' if \code{only_best} is TRUE, only the best motif would be returned (a data framw with a single row).
#'
#' @examples
#' res_screen <- screen_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#' head(res_screen)
#'
#' # only best match
#' screen_pwm(cluster_sequences_example, cluster_mat_example[, 1])
#'
#' # with r^2 metric
#' res_screen <- screen_pwm(sequences_example, response_mat_example[, 1], metric = "r2")
#' head(res_screen)
#'
#' @export
screen_pwm <- function(sequences, response, metric = NULL, dataset = all_motif_datasets(), motifs = NULL, parallel = getOption("prego.parallel", TRUE), only_best = FALSE, prior = 0.01, alternative = "two.sided", ...) {
    if (!is.null(motifs)) {
        dataset <- dataset %>% filter(motif %in% motifs)
    }

    if (!is.vector(response)) {
        if (is.matrix(response)) {
            cli_alert_info("The response is a matrix. The average will be used")
            response <- rowMeans(response, na.rm = TRUE)
        } else {
            cli_abort("The response must be a vector or a matrix")
        }
    }

    if (length(sequences) != length(response)) {
        cli_abort("The number of sequences and the number of rows in {.field response} do not match")
    }

    if (any(is.na(sequences))) {
        cli_abort("There are missing values in the sequences")
    }

    if (is.null(metric)) {
        metric <- ifelse(is_binary_response(response), "ks", "r2")
    }

    if (!is_binary_response(response) && metric == "ks") {
        cli_abort("The {.field metric} cannot be {.val ks} for a continuous response")
    }

    cli_alert_info("Performing PWM screening")
    res <- plyr::daply(dataset, "motif", function(x) {
        pwm <- compute_pwm(sequences, x, prior = prior, ...)
        if (metric == "ks") {
            return(suppressWarnings(ks.test(pwm[as.logical(response)], pwm[!as.logical(response)], alternative = alternative)$statistic))
        } else {
            return(cor(pwm, response)^2)
        }
        purrr::map_dbl(cluster_ids, function(cl) {
            suppressWarnings(ks.test(pwm[clusters == cl], pwm[clusters != cl], alternative = alternative)$statistic)
        })
    }, .parallel = parallel)

    res <- tibble::enframe(res, "motif", "score") %>%
        arrange(desc(score))

    if (only_best) {
        res <- res %>% slice(1)
    }

    return(res)
}
