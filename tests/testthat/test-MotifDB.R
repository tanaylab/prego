create_test_motif_db <- function(prior = 0.01) {
    # Create a simple motif database with two motifs
    motif_data <- tibble(
        motif = rep(c("motif1", "motif2"), each = 4),
        pos = rep(1, 8),
        A = c(0.7, 0.1, 0.2, 0.3, 0.1, 0.6, 0.2, 0.1),
        C = c(0.1, 0.7, 0.3, 0.2, 0.2, 0.1, 0.6, 0.2),
        G = c(0.1, 0.1, 0.3, 0.3, 0.6, 0.2, 0.1, 0.1),
        T = c(0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.6)
    )

    create_motif_db(motif_data, prior = prior)
}

# Test MotifDB class creation and validation
test_that("MotifDB object creation works with valid input", {
    motif_db <- create_test_motif_db()
    expect_s4_class(motif_db, "MotifDB")
    expect_true(validObject(motif_db))
})

test_that("MotifDB validates matrix dimensions", {
    motif_data <- tibble(
        motif = c("motif1"),
        pos = c(1),
        A = c(0.7),
        C = c(0.1),
        G = c(0.1),
        T = c(0.1)
    )
    expect_error(create_motif_db(motif_data), NA)
})

test_that("MotifDB validates prior constraints", {
    expect_error(create_test_motif_db(prior = 0))
    expect_error(create_test_motif_db(prior = 1))
    expect_error(create_test_motif_db(prior = -0.1))
})

test_that("MotifDB subsetting works correctly", {
    motif_db <- create_test_motif_db()

    # Test character subsetting
    subset_db <- motif_db["motif1"]
    expect_equal(ncol(subset_db@mat), 1)
    expect_equal(colnames(subset_db@mat), "motif1")
    expect_equal(names(subset_db@motif_lengths), "motif1")

    # Test numeric subsetting
    subset_db2 <- motif_db[1]
    expect_equal(subset_db@mat, subset_db2@mat)

    # Test multiple motif selection
    multi_subset <- motif_db[c("motif1", "motif2")]
    expect_equal(ncol(multi_subset@mat), 2)

    # Test error on invalid motif name
    expect_error(motif_db["invalid_motif"])
})

test_that("motif_db_to_dataframe conversion works correctly", {
    motif_db <- create_test_motif_db()
    df <- motif_db_to_dataframe(motif_db)

    # Check structure
    expect_equal(colnames(df), c("motif", "pos", "A", "C", "G", "T"))
    expect_equal(nrow(df), 8) # 2 motifs * 4 positions

    # Check probabilities sum to 1 (approximately due to floating point)
    sums <- df %>%
        rowwise() %>%
        mutate(sum = sum(c(A, C, G, T))) %>%
        pull(sum)
    expect_true(all(abs(sums - 1) < 1e-10))

    # Check that values are properly recovered
    # Test specific known values from create_test_motif_db
    first_pos <- df %>%
        filter(motif == "motif1", pos == 1) %>%
        select(A, C, G, T) %>%
        unlist()

    expect_equal(first_pos, c(A = 0.7, C = 0.1, G = 0.1, T = 0.1), tolerance = 1e-6)

    # Test round-trip conversion
    motif_db2 <- create_motif_db(df, prior = motif_db@prior)
    df2 <- motif_db_to_dataframe(motif_db2)

    # Compare all values in the dataframes
    expect_equal(
        df %>% arrange(motif, pos),
        df2 %>% arrange(motif, pos),
        tolerance = 1e-6
    )

    # Test with different priors
    df_prior_01 <- motif_db_to_dataframe(create_test_motif_db(prior = 0.01))
    df_prior_05 <- motif_db_to_dataframe(create_test_motif_db(prior = 0.05))

    # Values should be different with different priors
    expect_false(identical(
        df_prior_01 %>% arrange(motif, pos),
        df_prior_05 %>% arrange(motif, pos)
    ))
})

test_that("spatial factors validation works", {
    motif_db <- create_test_motif_db()

    # Test with valid spatial factors
    spat_factors <- matrix(1,
        nrow = 2,
        ncol = 3,
        dimnames = list(c("motif1", "motif2"), NULL)
    )
    expect_error(
        create_motif_db(
            motif_db_to_dataframe(motif_db),
            spat_factors = spat_factors
        ),
        NA
    )

    # Test with invalid dimensions
    invalid_spat_factors <- matrix(1, nrow = 3, ncol = 3)
    expect_error(
        create_motif_db(
            motif_db_to_dataframe(motif_db),
            spat_factors = invalid_spat_factors
        )
    )

    # Test with negative values
    invalid_values_spat_factors <- spat_factors
    invalid_values_spat_factors[1, 1] <- -1
    expect_error(
        create_motif_db(
            motif_db_to_dataframe(motif_db),
            spat_factors = invalid_values_spat_factors
        )
    )
})

test_that("reverse complement matrix is correctly computed", {
    motif_db <- create_test_motif_db()

    # For a single position, check if reverse complement is correct
    check_rc <- function(pos, nuc_idx) {
        forward <- motif_db@mat[pos * 4 - (4 - nuc_idx), ]
        rc_pos <- motif_db@motif_lengths[1] - pos + 1
        rc_nuc_idx <- c(4, 3, 2, 1)[nuc_idx] # A->T, C->G, G->C, T->A
        reverse <- motif_db@rc_mat[rc_pos * 4 - (4 - rc_nuc_idx), ]
        expect_equal(forward, reverse)
    }

    # Check all positions and nucleotides
    for (pos in 1:4) {
        for (nuc_idx in 1:4) {
            check_rc(pos, nuc_idx)
        }
    }
})

test_that("prior modification is equivalent to creating with new prior", {
    # Create original object
    motif_db <- create_test_motif_db(prior = 0.01)

    # Create new object with different prior
    new_prior <- 0.02
    motif_db_direct <- create_test_motif_db(prior = new_prior)

    # Modify original object's prior
    prior(motif_db) <- new_prior

    # Compare the two objects
    expect_equal(motif_db@mat, motif_db_direct@mat)
    expect_equal(motif_db@rc_mat, motif_db_direct@rc_mat)
    expect_equal(motif_db@motif_lengths, motif_db_direct@motif_lengths)
    expect_equal(motif_db@prior, motif_db_direct@prior)
    expect_equal(motif_db@spat_factors, motif_db_direct@spat_factors)
    expect_equal(motif_db@spat_bin_size, motif_db_direct@spat_bin_size)

    # Test with multiple prior changes
    priors <- c(0.01, 0.02, 0.05, 0.01)
    motif_db <- create_test_motif_db(prior = priors[1])

    for (new_prior in priors[-1]) {
        # Create new object directly
        motif_db_direct <- create_test_motif_db(prior = new_prior)

        # Modify existing object
        prior(motif_db) <- new_prior

        # Compare objects
        expect_equal(motif_db@mat, motif_db_direct@mat)
        expect_equal(motif_db@rc_mat, motif_db_direct@rc_mat)
        expect_equal(motif_db@prior, new_prior)

        # Check that probabilities still sum to 1 after conversion to dataframe
        df <- motif_db_to_dataframe(motif_db)
        sums <- df %>%
            rowwise() %>%
            mutate(sum = sum(c(A, C, G, T))) %>%
            pull(sum)
        expect_true(all(abs(sums - 1) < 1e-10))
    }
})

test_that("motif_lengths validation works", {
    motif_db <- create_test_motif_db()

    # Test that motif lengths match matrix dimensions
    expect_equal(
        length(motif_db@motif_lengths),
        ncol(motif_db@mat)
    )

    # Test that names match
    expect_equal(
        names(motif_db@motif_lengths),
        colnames(motif_db@mat)
    )

    # Test that lengths are positive
    expect_true(all(motif_db@motif_lengths > 0))
})

test_that("as.data.frame works correctly", {
    motif_db <- create_test_motif_db()

    # Test direct conversion
    df1 <- as.data.frame(motif_db)
    df2 <- as.data.frame(motif_db_to_dataframe(motif_db))
    # Results should be identical
    expect_equal(df1, df2)

    # Test with row.names and optional parameters (should ignore them)
    df3 <- as.data.frame(motif_db, row.names = letters[1:8], optional = TRUE)
    expect_equal(df1, df3)

    # Test structure
    expect_equal(colnames(df1), c("motif", "pos", "A", "C", "G", "T"))
    expect_equal(nrow(df1), 8) # 2 motifs * 4 positions
})
