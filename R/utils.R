#' Check if reponse is binary
#'
#' @param response a matrix of response values
#'
#' @noRd
is_binary_response <- function(response) {
    if (!is.matrix(response)) {
        response <- matrix(response, ncol = 1)
    }
    ncol(response) == 1 && all(response %in% c(0, 1))
}

#' Set parallel threads
#'
#' @description Set the number of parallel threads to use. prego uses the R function \code{doMC::registerDoMC} to register the parallelization.
#' By default, prego uses 80% of the number of available cores. The options are saved under 'prego.parallel' (should we use parallelization, logical) and 'prego.parallel.nc' (number of cores to use, integer).
#'
#' @param thread_num number of threads. use '1' for non parallel behavior
#'
#' @return None
#'
#' @examples
#' \donttest{
#' set_parallel(8)
#' }
#' @export
set_parallel <- function(thread_num = max(1, round(parallel::detectCores() * 0.8))) {
    if (thread_num <= 1) {
        options(prego.parallel = FALSE)
        cli_alert_info("Parallelization disabled.")
    } else {
        doMC::registerDoMC(thread_num)
        options(prego.parallel = TRUE)
        options(prego.parallel.nc = thread_num)
        cli_alert_info("Parallelization enabled. Using {.val {thread_num}} threads.")
    }
}
