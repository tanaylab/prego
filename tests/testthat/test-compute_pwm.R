res_samp <- regress_pwm.sample(cluster_sequences_example, cluster_mat_example[, 1], final_metric = "ks")

test_that("regression result is the same as the one from regress_pwm()$pred", {
    expect_equal(res_samp$pred, compute_pwm(cluster_sequences_example, pssm = res_samp$pssm, spat = res_samp$spat))
})
