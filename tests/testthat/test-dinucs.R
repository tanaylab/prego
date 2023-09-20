test_that("calc_sequences_dinuc_dist handles character vector input correctly", {
    sequences <- c("AACGT", "CGTAA", "GGCCA")

    # Expecting no error when valid sequences are passed
    expect_error(calc_sequences_dinuc_dist(sequences, size = 5), NA)

    # Expecting a data frame as output
    expect_true(is.data.frame(calc_sequences_dinuc_dist(sequences, size = 5)))
})

test_that("calc_sequences_dinuc_dist handles non-character input correctly", {
    non_char_input <- list("AACGT", "CGTAA", "GGCCA")

    # Expecting an error when non-character input is passed
    expect_error(calc_sequences_dinuc_dist(non_char_input, size = 5))
})

test_that("calc_sequences_dinuc_dist handles null size parameter correctly", {
    sequences <- c("AACGT", "CGTAA", "GGCCA")

    # When size is NULL, it should automatically take the max size from sequences
    res <- calc_sequences_dinuc_dist(sequences, size = NULL)
    expect_equal(nrow(res), 5)
})

test_that("calc_sequences_dinuc_dist handles sequences shorter than size correctly", {
    sequences <- c("AACGT", "CGTAA", "GG")

    # Expecting an error when any of the sequences is shorter than the specified size
    expect_error(calc_sequences_dinuc_dist(sequences, size = 5))
})

# Assuming that your Rcpp function is correctly computing the distributions, you might write tests like:
test_that("calc_sequences_dinuc_dist computes correct dinucleotide distribution", {
    sequences <- c("AA", "CC", "GG", "TT", "AC")

    # For simplicity, we're only testing a small set of sequences and a small size
    res <- calc_sequences_dinuc_dist(sequences, size = 2)

    # We expect 20% for each dinucleotide at position 1
    expect_equal(res$AA[1], 0.2)
    expect_equal(res$CC[1], 0.2)
    expect_equal(res$GG[1], 0.2)
    expect_equal(res$TT[1], 0.2)
    expect_equal(res$AC[1], 0.2)

    # At position 2, the dinucleotide frequencies should be NA since we can't form dinucleotides from a single base
    expect_true(is.na(res$AA[2]))
    expect_true(is.na(res$CC[2]))
    expect_true(is.na(res$GG[2]))
    expect_true(is.na(res$TT[2]))
    expect_true(is.na(res$AC[2]))
})
