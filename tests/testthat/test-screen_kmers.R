random_nuc <- function(n, len) purrr::map_chr(1:n, ~ paste(sample(c("C", "T", "G", "A"), size = len, replace = TRUE), collapse = ""))

get_seqs <- function(n, frac, pat, len = 50) {
    pat_seq <- paste0(random_nuc(round(n * frac), len), pat, random_nuc(round(n * frac), len))
    no_pat_seq <- random_nuc(round(n * (1 - frac)), nchar(pat_seq[1]))
    seqs <- c(pat_seq, no_pat_seq)
    return(seqs)
}

test_that("screen_kmers works with 1D", {
    withr::local_seed(60427)
    seqs <- get_seqs(1000, 0.5, "GATAAGA")
    resp <- c(rep(1, 500), rep(0, 500))

    res <- screen_kmers(seqs, resp, kmer_length = 7, min_cor = 0.08, min_n = 50, return_mat = TRUE)
    expect_true(res["GATAAGA", 1] > 0.5)
    expect_equal(sum(res[, 1] > 0.5), 1)
    expect_true(is.matrix(res))
    expect_true(!is.null(rownames(res)))
    expect_equal(ncol(res), 1)

    res_single <- screen_kmers(seqs, resp, kmer_length = 7, min_cor = 0.5, min_n = 50, return_mat = TRUE)
    expect_equal(nrow(res_single), 1)

    res_df <- screen_kmers(seqs, resp, kmer_length = 7, min_cor = 0.08, min_n = 50)
    expect_equal(colnames(res_df)[1:4], c("kmer", "max_r2", "avg_n", "avg_var"))
    expect_equal(unlist(res_df[, 5]), res[, 1], ignore_attr = TRUE)
    expect_true(res_df[res_df$kmer == "GATAAGA", "avg_n"] - 0.5 <= 1e-2)
})

test_that("screen_kmers works with 1D with gaps", {
    withr::local_seed(60427)
    seqs <- get_seqs(400, 0.5, "GATAAGA")
    resp <- c(rep(1, 200), rep(0, 200))

    res <- screen_kmers(seqs, resp, kmer_length = 7, min_cor = 0.08, min_n = 50, return_mat = TRUE, max_gap = 3)
    expect_true(res["GATAAGA", 1] > 0.5)
    expect_equal(sum(res[, 1] > 0.5), 1)
    expect_true(is.matrix(res))
    expect_true(!is.null(rownames(res)))
    expect_equal(ncol(res), 1)

    res_df <- screen_kmers(seqs, resp, kmer_length = 7, min_cor = 0.08, min_n = 50, max_gap = 3)
    expect_equal(colnames(res_df)[1:4], c("kmer", "max_r2", "avg_n", "avg_var"))
    expect_equal(unlist(res_df[, 5]), res[, 1], ignore_attr = TRUE)
    expect_true(res_df[res_df$kmer == "GATAAGA", "avg_n"] - 0.5 <= 1e-2)
})
