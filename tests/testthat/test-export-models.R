res_multi <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7,
    motif_num = 2
)

test_that("export multi motifs regression works", {
    export_fn <- tempfile()
    res_multi$export(export_fn)
    r <- load_multi_regression(export_fn)
    expect_equal(r$predict(cluster_sequences_example), res_multi$predict(cluster_sequences_example))
    expect_equal(r$predict_multi(cluster_sequences_example), res_multi$predict_multi(cluster_sequences_example))
    expect_equal(r$models[[1]]$predict(cluster_sequences_example), res_multi$models[[1]]$predict(cluster_sequences_example))

    expect_equal(r$models[[1]]$pssm, res_multi$models[[1]]$pssm)
    expect_equal(r$models[[1]]$spat, res_multi$models[[1]]$spat)
    expect_equal(r$model$pssm, res_multi$model$pssm)
    expect_equal(r$model$spat, res_multi$model$spat)
})
