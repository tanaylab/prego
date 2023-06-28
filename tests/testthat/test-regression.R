res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7, multi_kmers = FALSE
)

res_samp <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7,
    sample_for_kmers = TRUE
)

res_multi_kmers <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7, multi_kmers = TRUE
)


test_that("predict() and $pred are the same", {
    expect_equal(res$pred, res$predict(cluster_sequences_example))
})

test_that("predict() and $pred are the same when sample_for_kmers is TRUE", {
    expect_equal(res_samp$pred, res_samp$predict(cluster_sequences_example))
})

test_that("predict() and $pred are the same when multi_kmers is TRUE", {
    expect_equal(res_multi_kmers$pred, res_multi_kmers$predict(cluster_sequences_example))
})
