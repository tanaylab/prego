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
    expect_equal(dim(res), c(2, 3)) # 2 sequences and 3 kmers of length 2

    # Test 2: Check correct frequency calculation without gaps
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length)
    expect_equal(res[1, "AT"], 1, ignore_attr = TRUE) # 'AT' appears once in the first sequence
    expect_equal(res[2, "CG"], 1, ignore_attr = TRUE) # 'CG' appears once in the second sequence

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

    # Test 6: Check correct frequency calculation with max_gap = 1
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 4
    res <- kmer_matrix(sequences, kmer_length, max_gap = 1)
    expect_equal(res[1, "ATCG"], 1, ignore_attr = TRUE) # 'ATCG' appears once in the first sequence
    expect_equal(res[2, "NTCG"], 1, ignore_attr = TRUE) # 'NTCG' appears once in the second sequence
    expect_equal(res[2, "ANCG"], 1, ignore_attr = TRUE) # 'NTCG' appears once in the second sequence
    expect_equal(res[2, "ATNG"], 1, ignore_attr = TRUE) # 'ATNG' appears once in the second sequence
    expect_equal(res[2, "ATCN"], 1, ignore_attr = TRUE) # 'ATCN' appears once in the second sequence

    # Test 7: Check correct frequency calculation with max_gap = 2
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 4
    res <- kmer_matrix(sequences, kmer_length, max_gap = 2)
    expect_equal(res[1, "ATCG"], 1, ignore_attr = TRUE) # 'ATCG' appears once in the first sequence
    expect_equal(res[2, "NNCG"], 1, ignore_attr = TRUE) # 'NNTCG' appears once in the second sequence
    expect_equal(res[2, "ATNN"], 1, ignore_attr = TRUE) # 'ATNNG' appears once in the second sequence
    expect_equal(res[2, "ANNG"], 1, ignore_attr = TRUE) # 'ATCNN' appears once in the second sequence

    # Test 8: Check correct frequency calculation with max_gap = 0 (equivalent to no gap)
    sequences <- c("ATCG", "ATCG")
    kmer_length <- 2
    res <- kmer_matrix(sequences, kmer_length, max_gap = 0)
    expect_equal(res[1, "AT"], 1, ignore_attr = TRUE) # 'AT' appears once in the first sequence
    expect_equal(res[2, "CG"], 1, ignore_attr = TRUE) # 'CG' appears once in the second sequence
})

test_that("kmers_to_pssm handles single kmer correctly", {
    result <- kmers_to_pssm("ACGT", prior = 0.01)
    expect_equal(nrow(result), 4)
    expect_equal(ncol(result), 6)
    expect_equal(sum(result[result$pos == 1, c("A", "C", "G", "T")]), 1)
    expect_equal(sum(result[result$pos == 2, c("A", "C", "G", "T")]), 1)
    expect_equal(sum(result[result$pos == 3, c("A", "C", "G", "T")]), 1)
    expect_equal(sum(result[result$pos == 4, c("A", "C", "G", "T")]), 1)
})

test_that("kmers_to_pssm handles multiple kmers correctly", {
    result <- kmers_to_pssm(c("ACGT", "TGCA"), prior = 0.01)
    expect_equal(nrow(result), 8)
    expect_equal(ncol(result), 6)
    expect_equal(sum(result[result$pos == 1, c("A", "C", "G", "T")]), 2)
    expect_equal(sum(result[result$pos == 2, c("A", "C", "G", "T")]), 2)
    expect_equal(sum(result[result$pos == 3, c("A", "C", "G", "T")]), 2)
    expect_equal(sum(result[result$pos == 4, c("A", "C", "G", "T")]), 2)
})

test_that("kmers_to_pssm handles 'N' correctly", {
    result <- kmers_to_pssm("ACGN", prior = 0.01)
    expect_equal(nrow(result), 4)
    expect_equal(ncol(result), 6)
    expect_equal(sum(result[result$pos == 1, c("A", "C", "G", "T")]), 1)
    expect_equal(sum(result[result$pos == 2, c("A", "C", "G", "T")]), 1)
    expect_equal(sum(result[result$pos == 3, c("A", "C", "G", "T")]), 1)
    expect_equal(sum(result[result$pos == 4, c("A", "C", "G", "T")]), 1)
})
