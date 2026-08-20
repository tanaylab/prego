# local_serial_blas() runs on every PWM call, so it must never error. RhpcBLASctl
# reports NA - not NULL, and without raising - when the runtime it queries is
# absent, which is what omp_get_max_threads() does on macOS. Before this guard,
# that NA reached `NA > 1` and took down every calc_seq_pwm() call on that
# platform.

test_that("local_serial_blas survives RhpcBLASctl returning NA", {
    skip_if_not_installed("RhpcBLASctl")
    testthat::local_mocked_bindings(
        omp_get_max_threads = function(...) NA_integer_,
        blas_get_num_procs = function(...) NA_integer_,
        blas_set_num_threads = function(...) invisible(NULL),
        omp_set_num_threads = function(...) invisible(NULL),
        .package = "RhpcBLASctl"
    )
    expect_no_error(local_serial_blas())
})

test_that("calc_seq_pwm works when RhpcBLASctl reports NA thread counts", {
    skip_if_not_installed("RhpcBLASctl")
    testthat::local_mocked_bindings(
        omp_get_max_threads = function(...) NA_integer_,
        blas_get_num_procs = function(...) NA_integer_,
        blas_set_num_threads = function(...) invisible(NULL),
        omp_set_num_threads = function(...) invisible(NULL),
        .package = "RhpcBLASctl"
    )
    expect_no_error(calc_seq_pwm(c("ACGTACGT", "TGCATGCA"), MOTIF_DB))
})
