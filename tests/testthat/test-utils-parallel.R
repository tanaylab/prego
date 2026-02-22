test_that("set_parallel disables plyr parallel mode for one thread", {
    old_opts <- options()
    withr::defer(options(old_opts))

    set_parallel(1)

    expect_false(getOption("prego.parallel"))
    expect_equal(getOption("prego.parallel.nc"), 1)
})

test_that("set_parallel always records requested thread count", {
    old_opts <- options()
    withr::defer(options(old_opts))

    set_parallel(2)

    expect_equal(getOption("prego.parallel.nc"), 2)
})
