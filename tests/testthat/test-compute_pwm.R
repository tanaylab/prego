# logSumExp of a vector
log_sum_log_vec <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7,
)

s <- "CAGTAAAAGCTTTAATGCGTCTTGAGAGGGAGAGCATCAGCTTACAGAGCGAAGACCCCGAATGGCAAAACCCCGTCCCTTTTATGGAGAATTGCCCTCCGCCTCAGACACGTCGCTCCCTGATTGGCTGCAGCCCATCGGCCGAGTTGTCCTCACGGGGAAGGCAGAGCACATGGAGTGGAAAACTACCCCGGGCACATGCACAGATTACTTGTTTACTACTTAGAACACAGGATGTCAGCACCATCTTGTAATGGCGAATGTGAGGGCGGCTCCTCATACTTAGTTCCCTTTTTATGA"
pssm <- data.frame(pos = 0:14, A = c(
    0.240737825632095, 0.256109774112701,
    0.155907645821571, 0.001024218974635, 0.144046381115913, 0.844003200531006,
    0.562452733516693, 0.324804604053497, 0.248924136161804, 0.641762673854828,
    0.154276013374329, 0.0414229892194271, 0.275812089443207, 0.25,
    0.21894682943821
), C = c(
    0.230786353349686, 0.193948850035667,
    0.107831478118896, 0.12736551463604, 0.147198468446732, 0.00343735655769706,
    0.312835812568665, 0.190958619117737, 0.272343933582306, 0.00105959235224873,
    0.0164314024150372, 0.207740902900696, 0.181653186678886, 0.25,
    0.32295748591423
), G = c(
    0.248799309134483, 0.293831676244736,
    0.536404132843018, 0.001024218974635, 0.147198468446732, 0.152344271540642,
    0.0539040714502335, 0.00103273347485811, 0.304524749517441, 0.000942722195759416,
    0.346374750137329, 0.000539946369826794, 0.318346858024597, 0.25,
    0.229047849774361
), T = c(
    0.279676526784897, 0.256109774112701,
    0.199856758117676, 0.870585978031158, 0.561556696891785, 0.000215093168662861,
    0.0708073452115059, 0.483204007148743, 0.174207225441933, 0.356234937906265,
    0.482917785644531, 0.750296175479889, 0.224187895655632, 0.25,
    0.229047849774361
))

windows <- purrr::map_chr(1:(nchar(s) - 14), function(i) substr(s, i, i + 14))

test_that("regression predict() function works", {
    expect_equal(res$predict(cluster_sequences_example), res$pred)
})

test_that("compute_pwm works with logSumExp function", {
    r <- compute_pwm(s, pssm, func = "logSumExp")
    windows_r <- purrr::map_dbl(windows, compute_pwm, pssm, func = "logSumExp")

    expect_true(abs(log_sum_log_vec(purrr::map_dbl(windows, compute_pwm, pssm, func = "logSumExp")) - r) < 1e-4)
})


test_that("compute_pwm works with 'max' func", {
    r <- compute_pwm(s, pssm, func = "max")
    windows_r <- purrr::map_dbl(windows, compute_pwm, pssm, func = "max")
    expect_equal(max(windows_r), r)

    windows_log_sum_exp_r <- purrr::map_dbl(windows, compute_pwm, pssm, func = "logSumExp")
    expect_true(all(windows_r == windows_log_sum_exp_r))
})

test_that("compute_pwm handles N/* symmetrically on both strands", {
    # For a bidirectional motif, a window and its reverse-complement must score
    # identically. This held for plain ACGT but broke when the window contained
    # an N at an informative position: the forward strand used the column's
    # average log-prob while the reverse strand used a flat log(0.25).
    win <- "AAATAAAAAAAAAAA" # single motif-length window (nchar == nrow(pssm) == 15)
    win_n <- "AAANAAAAAAAAAAA" # N at position 4 (an informative, T-dominated column)
    # tolerance is well above float32 traversal noise (~1e-6) but far below the
    # ~1.6 nat gap the asymmetry produced.
    for (f in c("max", "logSumExp")) {
        expect_lt(abs(
            compute_pwm(win_n, pssm, func = f, bidirect = TRUE) -
                compute_pwm(rc(win_n), pssm, func = f, bidirect = TRUE)
        ), 1e-4)
        # sanity: the plain (N-free) window is strand-symmetric too
        expect_lt(abs(
            compute_pwm(win, pssm, func = f, bidirect = TRUE) -
                compute_pwm(rc(win), pssm, func = f, bidirect = TRUE)
        ), 1e-4)
    }
})

test_that("compute_pwm score of a sequence is independent of batch companions", {
    # A long sequence must score the same whether computed alone or in a batch
    # with a shorter sequence. Regression test for the scan window being capped
    # to nchar(sequences[1]) for every sequence.
    long_seq <- substr(s, 1, 100)
    short_seq <- substr(s, 1, 20)
    for (f in c("max", "logSumExp")) {
        alone <- compute_pwm(long_seq, pssm, func = f)
        with_short_first <- compute_pwm(c(short_seq, long_seq), pssm, func = f)[2]
        with_short_last <- compute_pwm(c(long_seq, short_seq), pssm, func = f)[1]
        expect_equal(with_short_first, alone, info = f)
        expect_equal(with_short_last, alone, info = f)
    }
})
