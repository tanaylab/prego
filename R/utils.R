#' Check if reponse is binary
#'
#' @param response a matrix of response values
#'
#' @noRd
is_binary_response <- function(response) {
    ncol(response) == 1 && all(response %in% c(0, 1))
}
