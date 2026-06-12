library(testthat)

make_pssm <- function(rows) {
    do.call(rbind, lapply(rows, function(r) {
        data.frame(A = r[1], C = r[2], G = r[3], T = r[4])
    }))
}

test_that("bind_pssms concatenates two PSSMs with no gap and no trimming", {
    p1 <- make_pssm(list(c(0.9, 0.05, 0.025, 0.025), c(0.1, 0.1, 0.7, 0.1)))
    p2 <- make_pssm(list(c(0.1, 0.7, 0.1, 0.1), c(0.025, 0.025, 0.05, 0.9)))

    out <- bind_pssms(p1, p2, gap = 0, bits_thresh = NULL)

    expect_equal(nrow(out), 4)
    expect_equal(colnames(out), c("A", "C", "G", "T", "pos"))
    expect_equal(out$pos, 0:3)
    expect_equal(out[1:2, c("A", "C", "G", "T")], p1, ignore_attr = TRUE)
    expect_equal(out[3:4, c("A", "C", "G", "T")], p2, ignore_attr = TRUE)
})

test_that("bind_pssms inserts a uniform spacer of the requested length", {
    p1 <- make_pssm(list(c(0.9, 0.05, 0.025, 0.025)))
    p2 <- make_pssm(list(c(0.025, 0.025, 0.05, 0.9)))

    out <- bind_pssms(p1, p2, gap = 3, bits_thresh = NULL)

    expect_equal(nrow(out), 5)
    # the three spacer rows are uniform 0.25
    spacer <- out[2:4, c("A", "C", "G", "T")]
    expect_true(all(as.matrix(spacer) == 0.25))
    # spacer contributes zero information
    expect_equal(bits_per_pos(spacer), rep(0, 3))
    expect_equal(out$pos, 0:4)
})

test_that("bind_pssms trims low information edges before binding", {
    # 4 positions, only the middle 2 are informative
    informative <- list(
        c(0.25, 0.25, 0.25, 0.25), # flat
        c(0.9,  0.05, 0.025, 0.025),
        c(0.025, 0.025, 0.05, 0.9),
        c(0.25, 0.25, 0.25, 0.25)  # flat
    )
    p <- make_pssm(informative)

    out <- bind_pssms(p, p, gap = 0, bits_thresh = 0.1)

    # 2 informative + 2 informative = 4 rows
    expect_equal(nrow(out), 4)
    expect_equal(out$pos, 0:3)
})

test_that("bind_pssms ignores extra columns like pos", {
    p1 <- make_pssm(list(c(0.9, 0.05, 0.025, 0.025)))
    p1$pos <- 0
    p2 <- make_pssm(list(c(0.025, 0.025, 0.05, 0.9)))
    p2$pos <- 0

    out <- bind_pssms(p1, p2, gap = 0, bits_thresh = NULL)

    expect_equal(nrow(out), 2)
    expect_equal(out$pos, 0:1)
})

test_that("bind_pssms accepts matrix input", {
    p1 <- matrix(c(0.9, 0.05, 0.025, 0.025), nrow = 1, dimnames = list(NULL, c("A", "C", "G", "T")))
    p2 <- matrix(c(0.025, 0.025, 0.05, 0.9), nrow = 1, dimnames = list(NULL, c("A", "C", "G", "T")))

    out <- bind_pssms(p1, p2, gap = 2, bits_thresh = NULL)

    expect_equal(nrow(out), 4)
    expect_true(all(as.matrix(out[2:3, c("A", "C", "G", "T")]) == 0.25))
})

test_that("bind_pssms validates inputs", {
    p1 <- make_pssm(list(c(0.9, 0.05, 0.025, 0.025)))

    expect_error(bind_pssms(p1, p1, gap = -1, bits_thresh = NULL))
    expect_error(bind_pssms(p1, p1, gap = 1.5, bits_thresh = NULL))
    expect_error(bind_pssms(p1, p1, gap = 0, bits_thresh = c(0.1, 0.2)))

    bad <- data.frame(A = 0.25, C = 0.25, G = 0.25, U = 0.25) # has U not T
    expect_error(bind_pssms(bad, p1, gap = 0, bits_thresh = NULL))
    expect_error(bind_pssms(p1, bad, gap = 0, bits_thresh = NULL))
})
