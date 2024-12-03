test_that("compute_pwm and calc_seq_pwm produce the same results", {
    test_seq <- "ACGTACGT"

    test_pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), nrow = 4, byrow = TRUE)
    colnames(test_pssm) <- c("A", "C", "G", "T")

    test_motif_db <- data.frame(
        motif = "test_motif",
        pos = 1:4,
        A = test_pssm[, 1],
        C = test_pssm[, 2],
        G = test_pssm[, 3],
        T = test_pssm[, 4]
    )


    test_mdb <- create_motif_db(test_motif_db)
    scores_new <- calc_seq_pwm(test_seq, test_mdb, bidirect = FALSE)

    # Compare with original implementation:
    scores_old <- compute_pwm(test_seq, test_pssm, bidirect = FALSE)
    expect_true(abs(scores_new - scores_old) < 1e-6)

    scores_new <- calc_seq_pwm(test_seq, test_mdb, bidirect = TRUE)

    # Compare with original implementation:
    scores_old <- compute_pwm(test_seq, test_pssm, bidirect = TRUE)
    expect_true(abs(scores_new - scores_old) < 1e-6)
})

test_that("compute_pwm and calc_seq_pwm produce the same results with spatial", {
    res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
        final_metric = "ks", spat_bin_size = 40,
        spat_num_bins = 7,
    )
    mdb1 <- create_motif_db(res$pssm %>% mutate(motif = "m1"), spat_factors = t(as.matrix(data.frame(m1 = res$spat$spat_factor))), spat_bin_size = 40)
    mdb2 <- create_motif_db(res$pssm %>% mutate(motif = "m1"))
    sq <- toupper(cluster_sequences_example)
    s <- calc_seq_pwm(sq, mdb1, bidirect = T)[, 1]
    s1 <- compute_pwm(sq, pssm = res$pssm, spat = res$spat, bidirect = T)
    sn <- calc_seq_pwm(sq, mdb2, bidirect = TRUE)[, 1]
    sn1 <- compute_pwm(sq, pssm = res$pssm)
    expect_true(all(abs(s - s1) < 1e-4))
    expect_true(all(abs(sn - sn1) < 1e-4))
})
