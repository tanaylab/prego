
# Function to count occurrences of a substring within a vector of strings
substr_count <- function(strings, substr) {
    sapply(strings, function(string) {
        stringr::str_count(string, substr)
    })
}

test_that("generate_kmers function works correctly", {
    # Test 1: Check if kmers of length 2 without gaps are correctly generated
    test_1 <- generate_kmers(2)
    expect_equal(length(test_1), 4^2)
    expect_true(all(nchar(test_1) == 2))
    expect_true(!any(grepl("N", test_1)))

    # Test 2: Check if kmers of length 3 with a single gap are correctly generated
    test_2 <- generate_kmers(3, min_gap = 1, max_gap = 1)
    expect_true(all(nchar(test_2) == 3))
    expect_true(all(substr_count(test_2, "N") == 1))

    # Test 3: Check if kmers of length 3 with a gap of 1 to 2 'N's are correctly generated
    test_3 <- generate_kmers(3, min_gap = 1, max_gap = 2)
    expect_true(all(nchar(test_3) == 3))
    expect_true(all(substr_count(test_3, "N") >= 1))
    expect_true(all(substr_count(test_3, "N") <= 2))

    # Test 4: Check if kmers of length 3 with a gap of 2 'N's only are correctly generated
    test_4 <- generate_kmers(3, min_gap = 2, max_gap = 2)
    expect_true(all(nchar(test_4) == 3))
    expect_true(all(substr_count(test_4, "N") == 2))
})

test_that("kmer_matrix function works correctly", {
    # Test 1: Check dimensions of the output matrix
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length)
    expect_equal(dim(res), c(2, 16)) # 2 sequences and 16 possible kmers of length 2

    # Test 2: Check correct frequency calculation without gaps
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length)
    expect_equal(res[1, "AT"], 1, ignore_attr = TRUE) # 'AT' appears once in the first sequence
    expect_equal(res[2, "CG"], 1, ignore_attr = TRUE) # 'CG' appears once in the second sequence

    # Test 3: Check correct frequency calculation with gaps
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length, min_gap = 1, max_gap = 1)
    expect_equal(res[1, "AN"], 1, ignore_attr = TRUE) # 'AN' appears twice in the first sequence considering N could be T/C/G
    expect_equal(res[2, "CN"], 1, ignore_attr = TRUE) # 'CN' appears once in the second sequence considering N could be G

    # Test 4: Check correct frequency calculation with from_range and to_range
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length, from_range = 1, to_range = 3)
    expect_equal(res[1, "AT"], 1, ignore_attr = TRUE) # 'AT' appears once in the first sequence
    expect_equal(res[2, "TC"], 1, ignore_attr = TRUE) # 'TC' appears once in the second sequence

    # Test 5: Check correct handling of rownames
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length, set_rownames = TRUE)
    expect_equal(rownames(res), sequences) # Row names should be the sequences
})
