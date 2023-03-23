res_multi <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7,
    motif_num = 2
)

test_that("predict() and $pred are the same for multi", {
    expect_equal(res_multi$pred, res_multi$predict(cluster_sequences_example))
})

test_that("predict_multi() works", {
    m_pred <- res_multi$predict_multi(cluster_sequences_example)
    expect_equal(res_multi$models[[1]]$predict(cluster_sequences_example), m_pred$e1)
})