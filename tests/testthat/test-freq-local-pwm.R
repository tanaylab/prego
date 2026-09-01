# Test fixtures: two motifs of different lengths, so padding and the
# per-motif NA tail are exercised everywhere.
freq_test_db <- function() {
    m4 <- data.frame(
        motif = "M4", pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.25),
        C = c(0.1, 0.7, 0.1, 0.25),
        G = c(0.1, 0.1, 0.7, 0.25),
        T = c(0.1, 0.1, 0.1, 0.25)
    )
    m6 <- data.frame(
        motif = "M6", pos = 1:6,
        A = c(0.9, 0.05, 0.4, 0.25, 0.1, 0.6),
        C = c(0.05, 0.9, 0.2, 0.25, 0.1, 0.2),
        G = c(0.03, 0.03, 0.2, 0.25, 0.7, 0.1),
        T = c(0.02, 0.02, 0.2, 0.25, 0.1, 0.1)
    )
    create_motif_db(rbind(m4, m6))
}

random_freqs <- function(m, seed = 42) {
    set.seed(seed)
    q <- matrix(stats::runif(m * 4, 0.05, 1), nrow = m, ncol = 4)
    q <- q / rowSums(q)
    colnames(q) <- c("A", "C", "G", "T")
    q
}

onehot_freqs <- function(seq) {
    q <- diag(4)[match(strsplit(seq, "")[[1]], c("A", "C", "G", "T")), ]
    colnames(q) <- c("A", "C", "G", "T")
    q
}

# Independent reference: the definition, written as nested loops.
brute_freq_local_pwm <- function(q, mdb, combine, bidirect) {
    n <- ncol(mdb@mat)
    m <- nrow(q)
    out <- matrix(NA_real_, n, m, dimnames = list(colnames(mdb@mat), NULL))
    one_strand <- function(mat, i, j, len) {
        s <- 0
        for (l in seq_len(len)) {
            p <- mat[((l - 1) * 4 + 1):(l * 4), i]
            if (combine == "sum") {
                s <- s + sum(q[j + l - 1, ] * p)
            } else {
                s <- s + log(sum(q[j + l - 1, ] * exp(p)))
            }
        }
        s
    }
    for (i in seq_len(n)) {
        len <- mdb@motif_lengths[i]
        for (j in seq_len(m - len + 1)) {
            fwd <- one_strand(mdb@mat, i, j, len)
            if (bidirect) {
                rev <- one_strand(mdb@rc_mat, i, j, len)
                out[i, j] <- log(exp(fwd) + exp(rev))
            } else {
                out[i, j] <- fwd
            }
        }
    }
    out
}

test_that("calc_freq_local_pwm on a one-hot matrix reproduces compute_local_pwm", {
    mdb <- freq_test_db()
    seq <- "ACGTACGTTGCAAGGTCCAT"
    q <- onehot_freqs(seq)
    db <- motif_db_to_dataframe(mdb)

    for (bidirect in c(FALSE, TRUE)) {
        res <- calc_freq_local_pwm(q, mdb, combine = "sum", bidirect = bidirect)
        for (motif in rownames(res)) {
            pssm <- db[db$motif == motif, c("A", "C", "G", "T")]
            ref <- compute_local_pwm(seq, pssm, bidirect = bidirect, prior = mdb@prior)
            valid <- seq_len(nchar(seq) - nrow(pssm) + 1)
            expect_equal(res[motif, valid], as.numeric(ref)[valid], tolerance = 1e-5)
        }
    }
})

test_that("both combine methods agree with compute_local_pwm when the input is one-hot", {
    mdb <- freq_test_db()
    seq <- "ACGTACGTTGCAAGGTCCAT"
    q <- onehot_freqs(seq)

    a <- calc_freq_local_pwm(q, mdb, combine = "sum", bidirect = TRUE)
    b <- calc_freq_local_pwm(q, mdb, combine = "multiply", bidirect = TRUE)
    expect_equal(a, b, tolerance = 1e-5)
})

test_that("calc_freq_local_pwm matches a brute-force loop on a random frequency matrix", {
    mdb <- freq_test_db()
    q <- random_freqs(30)

    for (combine in c("multiply", "sum")) {
        for (bidirect in c(FALSE, TRUE)) {
            res <- calc_freq_local_pwm(q, mdb, combine = combine, bidirect = bidirect)
            ref <- brute_freq_local_pwm(q, mdb, combine, bidirect)
            expect_equal(res, ref, tolerance = 1e-5)
        }
    }
})

test_that("a flat frequency matrix gives every motif the same floor under combine='multiply'", {
    mdb <- freq_test_db()
    q <- matrix(0.25, nrow = 20, ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T")))

    res <- calc_freq_local_pwm(q, mdb, combine = "multiply", bidirect = FALSE)
    for (motif in rownames(res)) {
        len <- mdb@motif_lengths[[motif]]
        expect_equal(
            unname(res[motif, 1]), len * log(0.25),
            tolerance = 1e-5
        )
    }

    # ...whereas the expected-log-likelihood floor is motif-dependent.
    res_sum <- calc_freq_local_pwm(q, mdb, combine = "sum", bidirect = FALSE)
    expect_false(isTRUE(all.equal(
        res_sum["M6", 1] / mdb@motif_lengths[["M6"]],
        res_sum["M4", 1] / mdb@motif_lengths[["M4"]]
    )))
})

test_that("positions where a motif does not fit are NA, per motif length", {
    mdb <- freq_test_db()
    m <- 12
    res <- calc_freq_local_pwm(random_freqs(m), mdb)

    expect_equal(dim(res), c(2L, m))
    expect_equal(rownames(res), c("M4", "M6"))
    for (motif in rownames(res)) {
        len <- mdb@motif_lengths[[motif]]
        expect_false(anyNA(res[motif, seq_len(m - len + 1)]))
        expect_true(all(is.na(res[motif, (m - len + 2):m])))
    }
})

test_that("the four input forms give identical results", {
    mdb <- freq_test_db()
    q1 <- random_freqs(15, seed = 1)
    q2 <- random_freqs(15, seed = 2)
    expected <- list(a = calc_freq_local_pwm(q1, mdb), b = calc_freq_local_pwm(q2, mdb))

    # list of m x 4
    expect_equal(calc_freq_local_pwm(list(a = q1, b = q2), mdb), expected)

    # list of 4 x m (transposed orientation)
    expect_equal(calc_freq_local_pwm(list(a = t(q1), b = t(q2)), mdb), expected)

    # 3-D array: n_mats x m x 4
    arr <- array(0, dim = c(2, 15, 4), dimnames = list(c("a", "b"), NULL, NULL))
    arr[1, , ] <- q1
    arr[2, , ] <- q2
    expect_equal(calc_freq_local_pwm(arr, mdb), expected)

    # stacked: n_mats x 4m, the calc_seq_pwm layout
    stacked <- rbind(a = as.vector(t(q1)), b = as.vector(t(q2)))
    expect_equal(calc_freq_local_pwm(stacked, mdb), expected)
})

test_that("calc_freq_local_pwm accepts a tidy motif data frame as well as a MotifDB", {
    mdb <- freq_test_db()
    q <- random_freqs(15)
    expect_equal(
        calc_freq_local_pwm(q, motif_db_to_dataframe(mdb)),
        calc_freq_local_pwm(q, mdb),
        tolerance = 1e-5
    )
})

test_that("calc_freq_local_pwm rejects frequencies that are not distributions", {
    mdb <- freq_test_db()
    q <- random_freqs(15)
    q[3, ] <- q[3, ] * 2
    expect_error(calc_freq_local_pwm(q, mdb), "sum to 1")

    expect_error(calc_freq_local_pwm(matrix(0.25, 5, 5), mdb), "4")
})

test_that("calc_freq_local_pwm errors when the motif is longer than the frequency matrix", {
    mdb <- freq_test_db()
    expect_error(calc_freq_local_pwm(random_freqs(3), mdb), "shorter")
})
