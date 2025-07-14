test_that("screen_local_pwm respects comparison operators", {
    # Simple sequences and corresponding PSSM (motif length = 4)
    sequences <- c("ACGTACGT", "TGCATGCA")

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1, # A
        0.1, 0.7, 0.1, 0.1, # C
        0.1, 0.1, 0.7, 0.1, # G
        0.1, 0.1, 0.1, 0.7 # T
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Compute local PWM scores (matrix)
    local_scores <- compute_local_pwm(
        sequences,
        pssm,
        bidirect = FALSE,
        prior = 0,
        return_list = FALSE
    )

    # Choose a threshold that will split the scores roughly in half
    thresh <- median(local_scores, na.rm = TRUE)

    # Expected positions for each operator
    expected_gt <- lapply(seq_len(nrow(local_scores)), function(i) {
        which(!is.na(local_scores[i, ]) & local_scores[i, ] > thresh)
    })
    expected_lt <- lapply(seq_len(nrow(local_scores)), function(i) {
        which(!is.na(local_scores[i, ]) & local_scores[i, ] < thresh)
    })

    # Results from screen_local_pwm
    res_gt <- screen_local_pwm(sequences, pssm, ">", thresh, bidirect = FALSE, prior = 0)
    res_lt <- screen_local_pwm(sequences, pssm, "<", thresh, bidirect = FALSE, prior = 0)

    expect_equal(res_gt, expected_gt)
    expect_equal(res_lt, expected_lt)

    # The two operator results should differ for the chosen threshold
    expect_false(identical(res_gt, res_lt))
})


test_that("screen_local_pwm throws an error for invalid operators", {
    seqs <- c("ACGTACGT")
    pssm <- matrix(
        c(
            0.25, 0.25, 0.25, 0.25,
            0.25, 0.25, 0.25, 0.25,
            0.25, 0.25, 0.25, 0.25,
            0.25, 0.25, 0.25, 0.25
        ),
        nrow = 4, byrow = TRUE
    )
    colnames(pssm) <- c("A", "C", "G", "T")

    expect_error(
        screen_local_pwm(seqs, pssm, "!=", 0, bidirect = FALSE, prior = 0),
        regexp = "must be one of"
    )
})


test_that("screen_local_pwm >= and <= behave as supersets of > and <", {
    sequences <- c("ACGTACGT", "TGCATGCA")
    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    local_scores <- compute_local_pwm(sequences, pssm, bidirect = FALSE, prior = 0, return_list = FALSE)
    thresh <- median(local_scores, na.rm = TRUE)

    res_gt <- screen_local_pwm(sequences, pssm, ">", thresh, bidirect = FALSE, prior = 0)
    res_ge <- screen_local_pwm(sequences, pssm, ">=", thresh, bidirect = FALSE, prior = 0)
    res_lt <- screen_local_pwm(sequences, pssm, "<", thresh, bidirect = FALSE, prior = 0)
    res_le <- screen_local_pwm(sequences, pssm, "<=", thresh, bidirect = FALSE, prior = 0)

    # >= should include all positions from >
    expect_true(all(vapply(seq_along(res_gt), function(i) all(res_gt[[i]] %in% res_ge[[i]]), logical(1))))

    # <= should include all positions from <
    expect_true(all(vapply(seq_along(res_lt), function(i) all(res_lt[[i]] %in% res_le[[i]]), logical(1))))

    # The union of < and > should be subset of <= and >= respectively (no gaps)
    expect_true(all(vapply(seq_along(res_gt), function(i) all(unique(c(res_gt[[i]], res_lt[[i]])) %in% unique(c(res_ge[[i]], res_le[[i]]))), logical(1))))
})

