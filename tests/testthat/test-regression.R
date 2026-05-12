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

test_that("return_all returns a list of candidate-kmer regressions sorted by val_score", {
    res_all <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
        final_metric = "ks", spat_bin_size = 40,
        spat_num_bins = 7, multi_kmers = TRUE,
        return_all = TRUE, max_cands = 3, match_with_db = FALSE
    )
    expect_type(res_all, "list")
    expect_gte(length(res_all), 1)
    # each element looks like a single-motif result
    expect_true(all(sapply(res_all, function(x) all(c("pssm", "spat", "pred", "seed_motif", "val_score", "predict") %in% names(x)))))
    # sorted by val_score descending
    val_scores <- sapply(res_all, function(x) x$val_score)
    expect_equal(val_scores, sort(val_scores, decreasing = TRUE))
    # named by seed_motif
    expect_equal(names(res_all), unname(sapply(res_all, function(x) x$seed_motif)))
    # predict() works on the full sequences for each element
    for (r in res_all) {
        expect_length(r$predict(cluster_sequences_example), length(cluster_sequences_example))
    }
})

test_that("return_all with sample_for_kmers refits each candidate on full data", {
    res_all <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
        final_metric = "ks", spat_bin_size = 40,
        spat_num_bins = 7, multi_kmers = TRUE,
        return_all = TRUE, max_cands = 3, match_with_db = FALSE,
        sample_for_kmers = TRUE, sample_frac = 0.3
    )
    expect_type(res_all, "list")
    expect_gte(length(res_all), 1)
    # Refitted on full data: pred should match predict() over the full sequences
    for (r in res_all) {
        expect_equal(r$pred, r$predict(cluster_sequences_example))
    }
})

test_that("return_all errors when incompatible with motif_num or motif", {
    expect_error(
        regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
            return_all = TRUE, motif_num = 2
        ),
        "return_all"
    )
    expect_error(
        regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
            return_all = TRUE, motif = "ACGTACGT"
        ),
        "return_all"
    )
    expect_error(
        regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
            return_all = TRUE, multi_kmers = FALSE
        ),
        "return_all"
    )
})
