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

# Pin the BLAS to a single thread for the remainder of the calling function.
#
# prego's PWM kernels (compute_pwm / compute_local_pwm / calc_seq_pwm) run a
# small BLAS `dgemm` per sequence INSIDE an RcppParallel/TBB `parallelFor` over
# sequences. With a multi-threaded BLAS (MKL, OpenBLAS), each TBB worker spawns
# its own BLAS thread team, so a single call opens `n_threads * n_blas_threads`
# OS threads - thousands on a many-core node. That oversubscribes the machine
# and, on a cluster job with a thread/process (cgroup pids) limit, makes the
# call FAIL with a thread-creation error. The per-sequence dgemm is tiny and the
# TBB loop over sequences is the right level of parallelism, so the inner BLAS
# must be serial. Restores the previous BLAS thread count on exit, so the rest
# of the session keeps its threaded BLAS. No-op if RhpcBLASctl is unavailable.
local_serial_blas <- function(.local_envir = parent.frame()) {
    if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        return(invisible(NULL))
    }
    # RhpcBLASctl reports NA - not NULL, and without erroring - when the runtime
    # it queries is absent; omp_get_max_threads() does exactly that on macOS. NA
    # has to be screened out here or the comparison errors, and this guard runs
    # on every PWM call, so an error in it takes the whole package down.
    worth_pinning <- function(n) !is.null(n) && !is.na(n) && n > 1
    old <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
    RhpcBLASctl::blas_set_num_threads(1)
    if (worth_pinning(old)) {
        withr::defer(RhpcBLASctl::blas_set_num_threads(old), envir = .local_envir)
    }

    # blas_set_num_threads() alone is a NO-OP on an OpenMP-linked OpenBLAS -
    # conda's `libopenblas-*-openmp` build, which is what the lab envs ship.
    # There openblas_set_num_threads() does nothing and omp_get_max_threads()
    # is what sizes the thread team, so pin that too or this guard silently
    # protects nothing. Measured on a 128-core node, extract_pwm() over
    # 3,000 x 500bp sequences and 20 motifs at set_parallel(16):
    #   unguarded          >600s (killed), 2049 threads
    #   blas pin only      >600s (killed), 2049 threads
    #   + omp pin (this)     4.8s,         1922 threads
    #   OMP_NUM_THREADS=1    0.47s,          17 threads
    # The omp pin is worth ~125x but cannot close the gap on its own: it writes
    # the CALLING thread's libgomp ICV, and the TBB workers are separate threads
    # that still inherit the global ICV. Only OMP_NUM_THREADS=1 in the
    # environment BEFORE R starts closes the rest, and no package can set that
    # for the user - libgomp reads it in an ELF constructor at process load, so
    # Sys.setenv() is always too late. Pinning it per worker thread from C++
    # would; see the NEWS entry for 0.0.11.
    old_omp <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NULL)
    if (worth_pinning(old_omp)) {
        RhpcBLASctl::omp_set_num_threads(1)
        withr::defer(RhpcBLASctl::omp_set_num_threads(old_omp), envir = .local_envir)
    }
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
