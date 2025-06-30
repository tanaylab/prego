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

# Test multi-variate response
test_that("regress_multiple_motifs works with multi-variate response", {
    # Test with 2 response variables from response_mat_example
    res_multi_resp <- regress_multiple_motifs(
        sequences_example[1:100],
        response_mat_example[1:100, 1:2],
        motif_num = 2,
        spat_bin_size = 40,
        spat_num_bins = 7,
        match_with_db = FALSE,
        verbose = FALSE
    )

    # Should complete without error and return valid structure
    expect_type(res_multi_resp, "list")
    expect_true("models" %in% names(res_multi_resp))
    expect_true("multi_stats" %in% names(res_multi_resp))
    expect_true("pred" %in% names(res_multi_resp))
    expect_equal(length(res_multi_resp$models), 2)

    # Predictions should work
    pred_test <- res_multi_resp$predict(sequences_example[1:10])
    expect_true(is.matrix(pred_test) || is.array(pred_test))
    expect_equal(nrow(pred_test), 10)
    expect_equal(ncol(pred_test), 2) # Should match number of response variables

    # Multi predictions should work
    pred_multi_test <- res_multi_resp$predict_multi(sequences_example[1:10])
    expect_s3_class(pred_multi_test, "data.frame")
    expect_equal(nrow(pred_multi_test), 10)
    expect_equal(ncol(pred_multi_test), 2) # Should match motif_num
    expect_true(all(c("e1", "e2") %in% colnames(pred_multi_test)))
})

test_that("regress_multiple_motifs works with single response variable", {
    # Test with single response variable (existing functionality)
    res_single_resp <- regress_multiple_motifs(
        sequences_example[1:100],
        response_mat_example[1:100, 1],
        motif_num = 2,
        spat_bin_size = 40,
        spat_num_bins = 7,
        match_with_db = FALSE,
        verbose = FALSE
    )

    # Should complete without error and return valid structure
    expect_type(res_single_resp, "list")
    expect_true("models" %in% names(res_single_resp))
    expect_true("multi_stats" %in% names(res_single_resp))
    expect_true("pred" %in% names(res_single_resp))
    expect_equal(length(res_single_resp$models), 2)

    # Predictions should work
    pred_test <- res_single_resp$predict(sequences_example[1:10])
    expect_type(pred_test, "double")
    expect_equal(length(pred_test), 10)

    # Multi predictions should work
    pred_multi_test <- res_single_resp$predict_multi(sequences_example[1:10])
    expect_s3_class(pred_multi_test, "data.frame")
    expect_equal(nrow(pred_multi_test), 10)
    expect_equal(ncol(pred_multi_test), 2) # Should match motif_num
    expect_true(all(c("e1", "e2") %in% colnames(pred_multi_test)))
})
