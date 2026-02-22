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

is_pkg_installed <- function(pkg) {
    nzchar(system.file(package = pkg))
}

#' Set parallel threads
#'
#' @description Set the number of parallel threads to use. prego uses
#' \code{doMC::registerDoMC} (when available) to register the
#' parallelization backend for \pkg{plyr}.
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
    options(prego.parallel.nc = thread_num)
    RcppParallel::setThreadOptions(numThreads = thread_num)

    if (thread_num <= 1) {
        options(prego.parallel = FALSE)
        return(invisible(NULL))
    }

    if (!requireNamespace("doMC", quietly = TRUE)) {
        cli_warn(c(
            "{.pkg doMC} is not installed.",
            "i" = "Falling back to non-parallel {.pkg plyr} execution.",
            "i" = "Install {.pkg doMC} on Unix-like systems to enable {.pkg plyr} parallelism."
        ))
        options(prego.parallel = FALSE)
        return(invisible(NULL))
    }

    doMC::registerDoMC(thread_num)
    options(prego.parallel = TRUE)
    invisible(NULL)
}

safe_llply <- function(.data, .fun, ..., .parallel = FALSE) {
    tryCatch(
        {
            plyr::llply(.data, .fun, ..., .parallel = .parallel)
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                # Retry without parallel execution
                plyr::llply(.data, .fun, ..., .parallel = FALSE)
            } else {
                # If it's a different error, re-throw it
                stop(e)
            }
        }
    )
}

safe_ldply <- function(.data, .fun, ..., .parallel = FALSE, .id = NA) {
    tryCatch(
        {
            plyr::ldply(.data, .fun, ..., .parallel = .parallel, .id = .id)
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                # Retry without parallel execution
                plyr::ldply(.data, .fun, ..., .parallel = FALSE, .id = .id)
            } else {
                # If it's a different error, re-throw it
                stop(e)
            }
        }
    )
}

safe_daply <- function(.data, .variables, .fun = NULL, ..., .progress = "none",
                       .inform = FALSE, .drop_i = TRUE, .drop_o = TRUE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::daply(.data, .variables, .fun, ...,
                .progress = .progress,
                .inform = .inform, .drop_i = .drop_i, .drop_o = .drop_o, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::daply(.data, .variables, .fun, ...,
                    .progress = .progress,
                    .inform = .inform, .drop_i = .drop_i, .drop_o = .drop_o, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_adply <- function(.data, .margins, .fun = NULL, ..., .expand = TRUE,
                       .progress = "none", .inform = FALSE, .parallel = FALSE, .id = NA) {
    tryCatch(
        {
            plyr::adply(.data, .margins, .fun, ...,
                .expand = .expand,
                .progress = .progress, .inform = .inform, .parallel = .parallel, .id = .id
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::adply(.data, .margins, .fun, ...,
                    .expand = .expand,
                    .progress = .progress, .inform = .inform, .parallel = FALSE, .id = .id
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_ddply <- function(.data, .variables, .fun = NULL, ..., .progress = "none",
                       .inform = FALSE, .drop = TRUE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::ddply(.data, .variables, .fun, ...,
                .progress = .progress,
                .inform = .inform, .drop = .drop, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::ddply(.data, .variables, .fun, ...,
                    .progress = .progress,
                    .inform = .inform, .drop = .drop, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_alply <- function(.data, .margins, .fun = NULL, ..., .expand = TRUE,
                       .progress = "none", .inform = FALSE, .parallel = FALSE, .dims = FALSE) {
    tryCatch(
        {
            plyr::alply(.data, .margins, .fun, ...,
                .expand = .expand,
                .progress = .progress, .inform = .inform, .parallel = .parallel, .dims = .dims
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::alply(.data, .margins, .fun, ...,
                    .expand = .expand,
                    .progress = .progress, .inform = .inform, .parallel = FALSE, .dims = .dims
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_dlply <- function(.data, .variables, .fun = NULL, ..., .progress = "none",
                       .inform = FALSE, .drop = TRUE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::dlply(.data, .variables, .fun, ...,
                .progress = .progress,
                .inform = .inform, .drop = .drop, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::dlply(.data, .variables, .fun, ...,
                    .progress = .progress,
                    .inform = .inform, .drop = .drop, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_aaply <- function(.data, .margins, .fun = NULL, ..., .expand = TRUE,
                       .progress = "none", .inform = FALSE, .drop = TRUE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::aaply(.data, .margins, .fun, ...,
                .expand = .expand,
                .progress = .progress, .inform = .inform, .drop = .drop, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::aaply(.data, .margins, .fun, ...,
                    .expand = .expand,
                    .progress = .progress, .inform = .inform, .drop = .drop, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_laply <- function(.data, .fun = NULL, ..., .progress = "none",
                       .inform = FALSE, .drop = TRUE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::laply(.data, .fun, ...,
                .progress = .progress,
                .inform = .inform, .drop = .drop, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::laply(.data, .fun, ...,
                    .progress = .progress,
                    .inform = .inform, .drop = .drop, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_maply <- function(.data, .fun = NULL, ..., .expand = TRUE,
                       .progress = "none", .inform = FALSE, .drop = TRUE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::maply(.data, .fun, ...,
                .expand = .expand,
                .progress = .progress, .inform = .inform, .drop = .drop, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::maply(.data, .fun, ...,
                    .expand = .expand,
                    .progress = .progress, .inform = .inform, .drop = .drop, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}

safe_mlply <- function(.data, .fun = NULL, ..., .expand = TRUE,
                       .progress = "none", .inform = FALSE, .parallel = FALSE) {
    tryCatch(
        {
            plyr::mlply(.data, .fun, ...,
                .expand = .expand,
                .progress = .progress, .inform = .inform, .parallel = .parallel
            )
        },
        error = function(e) {
            if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
                warning("Memory allocation error detected. Falling back to non-parallel execution...")
                plyr::mlply(.data, .fun, ...,
                    .expand = .expand,
                    .progress = .progress, .inform = .inform, .parallel = FALSE
                )
            } else {
                stop(e)
            }
        }
    )
}
