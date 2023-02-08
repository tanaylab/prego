test_that("Mask sequences by PWM", {
    seq <- c("ACCGGGGGTTTT", "TTTTTTTTTTTT", "TTTTTTTTACCG", "ACCGGGGGTTTT", "TTTTTTTTTTTT", "TTTTTTTTACCG")
    pwm <- data.frame(
        A = c(0.30, 0.10, 0.10, 0.30),
        C = c(0.10, 0.30, 0.30, 0.10),
        G = c(0.30, 0.30, 0.30, 0.30),
        T = c(0.30, 0.30, 0.30, 0.30)
    )
    res <- mask_sequences_by_pwm(seq, pwm, -4.6, pos_bits_thresh = 0)
    expect_equal(res, c("NNNNGNNNNNNN", "TTTTTTTTNNNN", "TTTTTTTTNNNN", "NNNNGNNNNNNN", "TTTTTTTTNNNN", "TTTTTTTTNNNN"))
})
