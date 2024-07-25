test_that("calc_sequences_dinucs works correctly", {
    sequences <- c("ATCG", "GCTA", "AATT")
    result <- calc_sequences_dinucs(sequences)

    expect_true("matrix" %in% class(result))
    expect_equal(dim(result), c(3, 16))
    expect_equal(colnames(result), c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"))

    expect_equal(result[1, ], c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0), ignore_attr = TRUE)
    expect_equal(result[2, ], c(0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0), ignore_attr = TRUE)
    expect_equal(result[3, ], c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1), ignore_attr = TRUE)

    expect_error(calc_sequences_dinucs(""))

    long_seq <- paste(rep("ATCG", 1000), collapse = "")
    result_long <- calc_sequences_dinucs(long_seq)
    expect_equal(dim(result_long), c(1, 16))
    expect_equal(sum(result_long), 3999) # Total dinucleotides in a sequence of length 4000

    identical_seqs <- rep("ATCG", 5)
    result_identical <- calc_sequences_dinucs(identical_seqs)
    expect_equal(dim(result_identical), c(5, 16))
    expect_true(all(t(result_identical) == result_identical[1, ]))

    lower_case <- "atcg"
    upper_case <- "ATCG"
    expect_equal(calc_sequences_dinucs(lower_case), calc_sequences_dinucs(upper_case))

    many_seqs <- replicate(10000, paste(sample(c("A", "T", "C", "G"), 10, replace = TRUE), collapse = ""))
    result_many <- calc_sequences_dinucs(many_seqs)
    expect_equal(dim(result_many), c(10000, 16))
})


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

test_that("calc_sequences_trinuc_dist handles character vector input correctly", {
    sequences <- c("AACGTT", "CGTAAG", "GGCCAA")

    # Expecting no error when valid sequences are passed
    expect_error(calc_sequences_trinuc_dist(sequences, size = 6), NA)

    # Expecting a data frame as output
    expect_true(is.data.frame(calc_sequences_trinuc_dist(sequences, size = 6)))
})

test_that("calc_sequences_trinuc_dist handles non-character input correctly", {
    non_char_input <- list("AACGTT", "CGTAAG", "GGCCAA")

    # Expecting an error when non-character input is passed
    expect_error(calc_sequences_trinuc_dist(non_char_input, size = 6))
})

test_that("calc_sequences_trinuc_dist handles null size parameter correctly", {
    sequences <- c("AACGTT", "CGTAAG", "GGCCAA")

    # When size is NULL, it should automatically take the max size from sequences
    res <- calc_sequences_trinuc_dist(sequences, size = NULL)
    expect_equal(nrow(res), 6)
})

test_that("calc_sequences_trinuc_dist handles sequences shorter than size correctly", {
    sequences <- c("AACGTT", "CGTAAG", "GG")

    # Expecting an error when any of the sequences is shorter than the specified size
    expect_error(calc_sequences_trinuc_dist(sequences, size = 6))
})

test_that("calc_sequences_trinuc_dist computes correct trinucleotide distribution", {
    sequences <- c("AAACCC", "GGGAAA", "TTTGGG", "CCCAAA", "AAAGGG")

    res <- calc_sequences_trinuc_dist(sequences, size = 6)

    expect_equal(res$AAA[1], 0.4)
    expect_equal(res$GGG[1], 0.2)
    expect_equal(res$TTT[1], 0.2)
    expect_equal(res$CCC[1], 0.2)

    expect_equal(res$AAC[2], 0.2)
    expect_equal(res$GGA[2], 0.2)
    expect_equal(res$TTG[2], 0.2)
    expect_equal(res$CCA[2], 0.2)
    expect_equal(res$AAG[2], 0.2)

    expect_true(all(is.na(res[5, -1])))
    expect_true(all(is.na(res[6, -1])))
})

test_that("calc_sequences_trinuc_dist handles edge cases correctly", {
    # Test with sequences exactly 3 nucleotides long
    short_sequences <- c("AAA", "CCC", "GGG", "TTT")
    res_short <- calc_sequences_trinuc_dist(short_sequences, size = 3)
    expect_equal(nrow(res_short), 3)
    expect_true(all(is.na(res_short[3, -1])))

    # Test with a single sequence
    single_sequence <- "AAAAAA"
    res_single <- calc_sequences_trinuc_dist(single_sequence, size = 6)
    expect_equal(res_single$AAA[1], 1)
    expect_true(all(res_single$AAA[2:4] == 1))
    expect_true(all(is.na(res_single[5:6, -1])))
})
