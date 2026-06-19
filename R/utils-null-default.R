#' Null-coalescing operator
#'
#' Returns `x` unless it is `NULL`, in which case it returns `y`. Defined
#' internally so that prego does not depend on the calling environment (or on
#' R >= 4.4, which exports `%||%` from base) for this operator. Without this,
#' uses of `%||%` inside functions dispatched to parallel workers (e.g. the
#' kmer-screen path in [regress_pwm()] / `learn_traj_prego()`) fail with
#' "could not find function \"%||%\"" on R 4.3 and earlier.
#'
#' @param x,y Values; `y` is returned only when `x` is `NULL`.
#' @return `x` if not `NULL`, otherwise `y`.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
