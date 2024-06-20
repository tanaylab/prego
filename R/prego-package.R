#' @keywords internal
"_PACKAGE"

#' @importFrom cli cli_abort cli_warn cli_alert cli_alert_info cli_ul cli_alert_success cli_h1 cli_h2 cli_h3 cli_progress_bar cli_progress_update cli_progress_done cli_progress_along cli_alert_warning
#' @importFrom tgstat tgs_cor tgs_dist tgs_matrix_tapply
#' @importFrom dplyr select mutate filter slice left_join right_join inner_join anti_join arrange desc pull group_by slice_max ungroup n_distinct distinct n everything any_of starts_with
#' @importFrom tibble tibble
#' @importFrom stats ks.test cor ecdf predict.lm lm predict approxfun quantile
#' @importFrom glue glue
#' @import ggplot2
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib prego
NULL
