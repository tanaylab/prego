library(testthat)

test_that("rc function works correctly", {
    # Test single sequence
    expect_equal(rc("ATCG"), "CGAT")

    # Test multiple sequences
    expect_equal(rc(c("ATCG", "GGCC", "TATA")), c("CGAT", "GGCC", "TATA"))

    # Test lowercase input
    expect_equal(rc("atcg"), "CGAT")

    # Test mixed case input
    expect_equal(rc("AtCg"), "CGAT")

    # Test empty string
    expect_equal(rc(""), "")

    # Test non-standard characters
    expect_equal(rc("AT-CG"), "CG-AT")

    # Test vector with empty string
    expect_equal(rc(c("ATCG", "", "GGCC")), c("CGAT", "", "GGCC"))

    # Test long sequence
    long_seq <- paste(rep("ATCG", 1000), collapse = "")
    expect_equal(rc(long_seq), paste(rep("CGAT", 1000), collapse = ""))

    # Test input with only non-standard characters
    expect_equal(rc("123"), "321")

    # Test NA input
    expect_true(is.na(rc(NA_character_)))

    # Test vector with NA
    expect_equal(rc(c("ATCG", NA, "GGCC")), c("CGAT", NA, "GGCC"))
})

test_that("rc function handles errors correctly", {
    # Test list input
    expect_error(rc(list("ATCG")), "The input should be a character vector")
})
