#' Export a motif regression model
#'
#' @param model a motif regression model, as returned by \code{regress_pwm} with \code{motif_num = 1}
#' @param fn a file name to save the model to
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
#'     final_metric = "ks", spat_bin_size = 40,
#'     spat_num_bins = 7
#' )
#' export_fn <- tempfile()
#' export_regression_model(export_fn)
#' r <- load_regression(export_fn)
#' }
#'
#' @export
export_regression_model <- function(model, fn) {
    r <- list(
        pssm = model$pssm,
        spat = model$spat,
        spat_min = model$spat_min,
        spat_max = model$spat_max,
        bidirect = model$bidirect,
        seq_length = model$seq_length
    )
    readr::write_rds(r, fn)
}

#' Load a motif regression model from a file
#'
#' @param fn file name
#'
#' @return a list with the following elements:
#'
#' \itemize{
#' \item{pssm:}{a data frame with the PSSM}
#' \item{spat:}{a data frame with the spatial profile}
#' \item{spat_min:}{a numeric value with the minimum value of the spatial profile}
#' \item{spat_max:}{a numeric value with the maximum value of the spatial profile}
#' \item{bidirect:}{a boolean value indicating whether the model is bidirectional}
#' \item{seq_length:}{a numeric value with the length of the sequences}
#' }
#'
#' @examples
#' \dontrun{
#' res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
#'     final_metric = "ks", spat_bin_size = 40,
#'     spat_num_bins = 7
#' )
#' export_fn <- tempfile()
#' export_regression_model(export_fn)
#' r <- load_regression(export_fn)
#' }
#'
#' @export
load_regression_model <- function(fn) {
    r <- readr::read_rds(fn)
    r$predict <- function(x, ...) {
        compute_pwm(x, r$pssm, spat = r$spat, bidirect = r$bidirect, spat_min = r$spat_min, spat_max = r$spat_max - 1, ...)
    }
    return(r)
}

#' Export a multiple motif regression model
#'
#' @param reg a multiple motif regression model, as returned by \code{regress_pwm} with \code{motif_num > 1}
#' @param fn a file name to save the model to
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' res_multi <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
#'     final_metric = "ks", spat_bin_size = 40,
#'     spat_num_bins = 7,
#'     motif_num = 2
#' )
#' export_fn <- tempfile()
#' export_multi_regression(res_multi, export_fn)
#'
#' # loading can be done by:
#' r <- load_multi_regression(export_fn)
#' }
#'
#' @export
export_multi_regression <- function(reg, fn) {
    export_model <- function(pssm, spat, spat_min, spat_max, bidirect, seq_length) {
        list(
            pssm = pssm,
            spat = spat,
            spat_min = spat_min,
            spat_max = spat_max,
            bidirect = bidirect,
            seq_length = seq_length
        )
    }

    models <- purrr::map(reg$models, ~ export_model(.x$pssm, .x$spat, .x$spat_min, .x$spat_max, .x$bidirect, .x$seq_length))

    model <- reg$model

    new_reg <- list(
        models = models,
        model = model,
        spat_min = reg$spat_min,
        spat_max = reg$spat_max,
        bidirect = reg$bidirect,
        spat_bin_size = reg$spat_bin_size,
        seq_length = reg$seq_length,
        motif_num = length(models)
    )

    readr::write_rds(new_reg, fn)
}

#' Load a multiple motif regression model from a file
#'
#' @param fn file name
#'
#' @return a list with the following elements:
#'
#' \itemize{
#'  \item{models: }{a list of models.}
#'  \item{model: }{the combined model.}
#'  \item{spat_min: }{the minimum spatial position.}
#'  \item{spat_max: }{the maximum spatial position.}
#'  \item{bidirect: }{whether the model is bidirectional.}
#'  \item{spat_bin_size: }{the spatial bin size.}
#'  \item{seq_length: }{the sequence length.}
#'  \item{motif_num: }{the number of motifs.}
#'  \item{predict: }{a function to predict the response.}
#'  \item{predict_multi: }{a function to predict the response for each motif.}
#' }
#'
#' @examples
#' \dontrun{
#' res_multi <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
#'     final_metric = "ks", spat_bin_size = 40,
#'     spat_num_bins = 7,
#'     motif_num = 2
#' )
#' tmp <- tempfile()
#' res_multi$export(tmp)
#' r <- load_multi_regression(tmp)
#' }
#'
#' @export
load_multi_regression <- function(fn) {
    r <- readr::read_rds(fn)
    r$models <- purrr::map(r$models, function(.x) {
        .x$predict <- function(x, ...) {
            compute_pwm(x, .x$pssm, spat = .x$spat, bidirect = .x$bidirect, spat_min = .x$spat_min, spat_max = .x$spat_max - 1, ...)
        }
        return(.x)
    })
    r$predict <- function(x, ...) {
        e <- lapply(1:r$motif_num, function(i) r$models[[i]]$predict(x, ...))
        e <- as.data.frame(e)
        colnames(e) <- paste0("e", 1:r$motif_num)
        predict(r$model, e, ...)
    }
    r$predict_multi <- function(x, parallel = getOption("prego.parallel", FALSE)) {
        e <- plyr::llply(r$models, function(.x) .x$predict(x), .parallel = parallel) %>%
            do.call(cbind, .) %>%
            as.data.frame()
        colnames(e) <- paste0("e", seq_along(r$models))
        return(e)
    }
    return(r)
}
